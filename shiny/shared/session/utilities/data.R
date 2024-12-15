# functions to load, parse, and filter bin data
protAtacLoadCreate <- "asNeeded" # for convenience when debugging load functions

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

    # for memory management, ensure that only one sourceId is loaded at a time
    workingSourceId <- isolate(binsWorkingSourceId())
    if(!is.null(workingSourceId) && workingSourceId != sourceId){
        binsCache$clear()
        binsWorkingSourceId(sourceId)
    } else if(is.null(workingSourceId)){
        binsWorkingSourceId(sourceId)
    }

    # load the analysis data into cache as needed
    binsCache$get(
        "paBinData",
        keyObject = list(
            sourceId = sourceId
        ),
        permanent = FALSE, # don't resave to disk, is large, complex object, slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){

            # load the bin data from the pipeline data package
            startSpinner(session, message = "loading bin counts.")
            bd <- readRDS(getSourceFilePath(sourceId, "binCounts"))

            # create a shared bin identifier for visualization merge operations
            startSpinner(session, message = "loading bin counts..")
            for(refType in refTypes) bd$bins[[refType]]$binI <- 1:nrow(bd$bins[[refType]])
            stopSpinner(session)
            bd
        }
    )$value
}

# load and format insert size data (from atac/collate action)
paInsertSizes <- function(sourceId){

    # load the analysis data into cache as needed
    insertSizesCache$get(
        "paInsertSizes",
        keyObject = list(
            sourceId = sourceId
        ),
        permanent = FALSE, # don't resave to disk, is large, complex object, slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){
            startSpinner(session, message = "loading inserts")
            isd <- readRDS(getSourceFilePath(sourceId, "insertSizes"))
            stopSpinner(session)
            # TODO: additional post-processing?
            isd
        }
    )$value
}

# load and format genome bins and read counts (from atat/collate action)
paScores <- function(sourceId){

    # for memory management, ensure that only one sourceId is loaded at a time
    workingSourceId <- isolate(scoresWorkingSourceId())
    if(!is.null(workingSourceId) && workingSourceId != sourceId){
        binsCache$clear()
        scoresWorkingSourceId(sourceId)
    } else if(is.null(workingSourceId)){
        scoresWorkingSourceId(sourceId)
    }

    # load the analysis data into cache as needed
    scoresCache$get(
        "paScores",
        keyObject = list(
            sourceId = sourceId
        ),
        permanent = FALSE, # don't resave to disk, is large, complex object, slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){
            startSpinner(session, message = "loading scores")
            sd <- readRDS(getSourceFilePath(sourceId, "scores"))
            stopSpinner(session)
            sd$reverseStageTypes <- {
                reversed <- list()
                for (name in names(sd$stageTypes)) {
                    for (value in sd$stageTypes[[name]]) {
                        reversed[[value]] <- name
                    }
                }
                reversed
            }
            # TODO: additional post-processing?
            sd
        }
    )$value
}

# load and format segmentation data (from atac/segment action)
paSegmentation <- function(sourceId){

    # load the analysis data into cache as needed
    segmentationCache$get(
        "paSegmentation",
        keyObject = list(
            sourceId = sourceId
        ),
        permanent = FALSE, # don't resave to disk, is large, complex object, slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){
            startSpinner(session, message = "loading segmentation")
            seg <- readRDS(getSourceFilePath(sourceId, "segmentation"))
            # TODO: additional post-processing?
            seg
        }
    )$value
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
        # create = protAtacLoadCreate,
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
            # for(key in names(sampleScores)) sampleScores[[key]]$score <- NULL
            # for(key in names(aggregateScores$by_stage)) aggregateScores$by_stage[[key]]$score <- NULL
            # for(key in names(aggregateScores$by_stageType)) aggregateScores$by_stageType[[key]]$score <- NULL
            # aggregateScores$stageType_delta$score <- NULL 
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
