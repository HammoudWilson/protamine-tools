# functions to load, parse, and filter bin data

# handle bin exclusions, including at gc extremes
isExcludedBin <- function(excluded, pct_gc){
    excluded == 1 | pct_gc < gcLimits[1] | pct_gc > gcLimits[2]
}
isIncludedAutosomeBin <- function(excluded, pct_gc, nAlleles){
    !isExcludedBin(excluded, pct_gc) & nAlleles == 2
}
getIncudedAutosomeBins <- function(bins){
    isIncludedAutosomeBin(bins$excluded, bins$pct_gc, bins$nAlleles)
}

# load and format genome bins and read counts (from atat/collate action)
paBinData <- function(sourceId){
    startSpinner(session, message = "loading bin counts")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "binCounts", 
        ttl = CONSTANTS$ttl$month, 
        postProcess = function(bd){
            for(refType in refTypes) bd$bins[[refType]]$binI <- 1:nrow(bd$bins[[refType]])
            bd
        }
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# load and format insert size data (from atac/collate action)
paInsertSizes <- function(sourceId){
    startSpinner(session, message = "loading insert sizes")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "insertSizes", 
        ttl = CONSTANTS$ttl$month, 
        postProcess = NULL
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# load and format genome bins and read counts (from atat/collate action)
paScores <- function(sourceId){
    startSpinner(session, message = "loading scores")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "scores", 
        ttl = CONSTANTS$ttl$month, 
        postProcess = function(sd){
            startSpinner(session, message = "scores post-processing")
            sd$reverseStageTypes <- {
                reversed <- list()
                for (name in names(sd$stageTypes)) {
                    for (value in sd$stageTypes[[name]]) {
                        reversed[[value]] <- name
                    }
                }
                reversed
            }
            sd
        }
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# load and format TSS data
paTssFragsData <- function(sourceId){
    startSpinner(session, message = "loading TSS inserts")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "tssFrags", 
        ttl = CONSTANTS$ttl$month
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# GC Residual Z-Score analysis, depends on pipeline bin data and user-selected GC bias models
#----------------------------------------------------------------------
# analyze and aggregate distributions of different bin scores
# all scores are expected to be one a comparable scale between samples, including
replaceNaN <- function(x){
    x[is.nan(x)] <- NA
    x
}
analyzeScoreDist <- function(bins, gcLimits, scores, scoreType){
    if(scoreType$log10) scores <- log10(pmax(scoreType$minValue, scores)) # prevent log(0) and impossible values
    scores_wrk <- scores[getIncudedAutosomeBins(bins)] # only use good autosomal bins to analyze score distributions
    dist <- data.table(x = as.integer(round(scores_wrk / scoreType$distUnit)) * scoreType$distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / replaceNaN(sum(y, na.rm = TRUE))]
    mu <- replaceNaN(mean(scores_wrk, na.rm = TRUE))
    sd <- replaceNaN(sd(  scores_wrk, na.rm = TRUE))
    median <- replaceNaN(median(scores_wrk, na.rm = TRUE))
    list(
        score    = scores, # score are dropped downstream when not included in scoreType
        dist     = dist,
        median   = median,
        mean     = mu,
        sd       = sd,
        z        = if("z" %in% scoreType$include) (scores - median) / sd else NULL,         # for normal/parametric scores
        quantile = if("quantile" %in% scoreType$include) ecdf(scores_wrk)(scores) else NULL # for non-parametric scores
    )
}
analyzeSampleScores <- function(bd, gcLimits, scoreFn, scoreType, ...){
    x <- lapply(bd$samples$sample_name, function(sample_name){
        startSpinner(session, message = paste("loading gcrz", sample_name))
        analyzeScoreDist(bd$bins$genome, gcLimits, scoreFn(bd, sample_name, ...), scoreType)
    })
    names(x) <- bd$samples$sample_name
    x
}
aggregateAndAnalyzeScores <- function(bins, gcLimits, sampleScores, sample_names, scoreType){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        sampleScores[[sample_name]]$score
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(bins, gcLimits, x[, replaceNaN(rowMeans(.SD, na.rm = TRUE))], scoreType)
}
aggregateSampleScores <- function(bd, gcLimits, stageTypes, sampleScores, scoreType){

    # aggregate scores by spermatid stage
    allStages <- unique(bd$samples$stage)
    by_stage <- lapply(allStages, function(stage_){
        startSpinner(session, message = paste("loading gcrz", stage_))
        sample_names <- bd$samples[stage == stage_, sample_name]
        aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, scoreType)
    })
    names(by_stage) <- allStages

    # aggregate scores by spermatid stage type (round vs. elong)
    by_stageType <- lapply(names(stageTypes), function(stageType){
        startSpinner(session, message = paste("loading gcrz", stageType))
        sample_names <- bd$samples[stage %in% stageTypes[[stageType]], sample_name]
        aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, scoreType)
    })
    names(by_stageType) <- names(stageTypes)

    # calculate the difference in scores between round and elongated spermatids to assess spermiogenesis trajectory
    startSpinner(session, message = paste("loading gcrz", "stageType_delta"))
    stageType1 <- names(stageTypes)[1]
    stageType2 <- names(stageTypes)[2]
    stageType_delta <- analyzeScoreDist(
        bd$bins$genome, 
        gcLimits, 
        by_stageType[[stageType1]]$score - by_stageType[[stageType2]]$score, # typically round - elong to give positive values in early stages
        scoreType
    )
    list(
        by_stage        = by_stage,
        by_stageType    = by_stageType,
        stageType_delta = stageType_delta
    )
}
gcResidualsZ <- function(sourceId, gcBiasModels = NULL){
    if(is.null(gcBiasModels)) gcBiasModels <- app$normalizeGC$getGcBiasModels_externalCall(sourceId)
    gcrzCache$get(
        "gcResidualsZ",
        keyObject = list(
            sourceId = sourceId,
            gcBiasModels = gcBiasModels
        ),
        permanent = TRUE,
        from = "ram",
        create = "asNeeded",
        createFn = function(...){
            startSpinner(session, message = "clearing gcrz")
            gcrzCache$clear()
            startSpinner(session, message = "loading gcrz")
            bd <- paBinData(sourceId)
            sd <- paScores(sourceId)
            scoreTypeName <- "gcrz"
            scoreType <- scoreTypes$sample[[scoreTypeName]]
            sampleScores <- analyzeSampleScores(bd, gcLimits, function(bd, sample_name){
                binCounts <- bd$binCounts$genome[, "all_inserts", sample_name]
                x <- app$normalizeGC$getBinZScore(sourceId, sample_name, binCounts, bd$bins$genome$pct_gc, bd$bins$genome$nAlleles)
                x[abs(x) == Inf] <- NA
                x
            }, scoreType)
            aggregateScores <- aggregateSampleScores(bd, sd$gcLimits, sd$stageTypes, sampleScores, scoreType)
            x <- list(
                sampleScores    = sampleScores,
                aggregateScores = aggregateScores
            )
            stopSpinner(session)
            x
        }
    )$value
}

# functions to load, parse, and filter bin data
getScoreLevel <- function(scoreTypeName){
         if(scoreTypeName %in% names(scoreTypes$genome)) 'genome'
    else if(scoreTypeName %in% names(scoreTypes$sample)) 'sample'
    else "NA"
}
getScoreType <- function(sourceId, scoreTypeName){
    scoreLevel <- getScoreLevel(scoreTypeName)
    scoreTypes[[scoreLevel]][[scoreTypeName]]
}
getTypedStages <- function(sourceId) unlist(paScores(sourceId)$stageTypes)
getStageTypesByStage <- function(sourceId, stages) unlist(paScores(sourceId)$reverseStageTypes[stages])

#----------------------------------------------------------------------
# sample-level score structure summary and associated score object retrieval
#----------------------------------------------------------------------
getGenomeScores <- function(sourceId, scoreTypeName){ # returns a single genome-level score object
    x <- list(paScores(sourceId)$scores$genome[[scoreTypeName]])
    names(x) <- scoreTypeName
    x
}
getSampleScoresList <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects based on GC normalization
    if(scoreTypeName == "gcrz"){
        gcResidualsZ(sourceId)
    } else {
        paScores(sourceId)$scores$sample[[scoreTypeName]]
    }
}
getSampleScores <- function(sourceId, scoreTypeName, samples){ # returns a list of sample-level score objects
    getSampleScoresList(sourceId, scoreTypeName)$sampleScores[samples$sample_name]
}
getStageScores <- function(sourceId, scoreTypeName, samples){ # returns a list of stage-level score objects matching a list of samples
    getSampleScoresList(sourceId, scoreTypeName)$aggregateScores$by_stage[unique(samples$stage)]
}
getStageTypeScores <- function(sourceId, scoreTypeName, samples){ # returns a list of stageType-level score objects matching a list of samples
    stageTypes <- getStageTypesByStage(sourceId, samples$stage)
    getSampleScoresList(sourceId, scoreTypeName)$aggregateScores$by_stageType[unique(stageTypes)]
}
getStageTypeDeltaScores <- function(sourceId, scoreTypeName, cleanDist = FALSE){ # returns a single score object for the delta between stage types (or a genome-level score object)
    list(
        stageType_delta = if(getScoreLevel(scoreTypeName) == "sample") {
            getSampleScoresList(sourceId, scoreTypeName)$aggregateScores$stageType_delta
        } else {
            x <- paScores(sourceId)$scores$genome[[scoreTypeName]]
            if(scoreTypeName == "txn" && cleanDist) x$dist <- x$dist[x %between% scoreTypes$genome$txn$valueLim]
            x
        }
    )
}
getSeriesAggScores <- function(sourceId, scoreTypeName, samples, config){
    switch(
        config$Aggregate_By,
        sample      = getSampleScores(sourceId, scoreTypeName, samples),
        stage       = getStageScores(sourceId, scoreTypeName, samples),
        stage_type  = getStageTypeScores(sourceId, scoreTypeName, samples)
    )
}
