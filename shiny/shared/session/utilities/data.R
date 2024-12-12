# functions to load, parse, and filter bin data
protAtacLoadCreate <- "once" # for convenience when debugging load functions

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

# load and format insert size data (from atat/collate action)
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

            # TODO: remove this and constants.R object once pipeline is re-run?
            sd$scoreTypes <- scoreTypes

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

# analyze and aggregate distributions of different bin scores
# all scores are expected to be one a comparable scale between samples, including
replaceNaN <- function(x){
    x[is.nan(x)] <- NA
    x
}
analyzeScoreDist <- function(bins, gcLimits, scores, distUnit = 0.1){
    scores_wrk <- scores[getIncudedAutosomeBins(bins)] # only use good autosomal bins to analyze score distributions
    dist <- data.table(x = as.integer(round(scores_wrk / distUnit)) * distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / replaceNaN(sum(y, na.rm = TRUE))]
    mu <- replaceNaN(mean(scores_wrk, na.rm = TRUE))
    sd <- replaceNaN(sd(  scores_wrk, na.rm = TRUE))
    list(
        score    = scores,
        dist     = dist,
        mean     = mu,
        sd       = sd,
        z        = (scores - mu) / sd
        # ,      # for normal/parametric scores
        # quantile = ecdf(scores_wrk)(scores) # for non-parametric scores
    )
}
analyzeSampleScores <- function(bd, gcLimits, scoreFn, distUnit = 0.1, ...){
    x <- lapply(bd$samples$sample_name, function(sample_name){
        startSpinner(session, message = paste("loading gcrz", sample_name))
        analyzeScoreDist(bd$bins$genome, gcLimits, scoreFn(bd, sample_name, ...), distUnit)
    })
    names(x) <- bd$samples$sample_name
    x
}
aggregateAndAnalyzeScores <- function(bins, gcLimits, sampleScores, sample_names, distUnit = 0.1){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        sampleScores[[sample_name]]$score
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(bins, gcLimits, x[, replaceNaN(rowMeans(.SD, na.rm = TRUE))], distUnit)
}
aggregateSampleScores <- function(bd, gcLimits, stageTypes, sampleScores, distUnit = 0.1){

    # aggregate scores by spermatid stage
    allStages <- unique(bd$samples$stage)
    by_stage <- lapply(allStages, function(stage_){
        startSpinner(session, message = paste("loading gcrz", stage_))
        sample_names <- bd$samples[stage == stage_, sample_name]
        aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, distUnit)
    })
    names(by_stage) <- allStages

    # aggregate scores by spermatid stage type (round vs. elong)
    by_stageType <- lapply(names(stageTypes), function(stageType){
        startSpinner(session, message = paste("loading gcrz", stageType))
        sample_names <- bd$samples[stage %in% stageTypes[[stageType]], sample_name]
        aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, distUnit)
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
        distUnit
    )
    list(
        by_stage        = by_stage,
        by_stageType    = by_stageType,
        stageType_delta = stageType_delta
    )
}
gcResidualsZ <- function(sourceId, gcBiasModels = NULL){
    if(is.null(gcBiasModels)) app$normalizeGC$getGcBiasModels_externalCall(sourceId)
    gcrzCache$get(
        "gcResidualsZ",
        keyObject = list(
            sourceId = sourceId,
            gcBiasModels = gcBiasModels
        ),
        permanent = FALSE, # don't resave to disk, is large, complex object, slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){
            startSpinner(session, message = "loading gcrz")
            bd <- paBinData(sourceId)
            sd <- paScores(sourceId)
            scoreTypeName <- "gcrz"
            scoreType <- sd$scoreTypes$sample[[scoreTypeName]]
            sampleScores <- analyzeSampleScores(bd, gcLimits, function(bd, sample_name){
                binCounts <- rowSums(bd$binCounts$genome[, , sample_name], na.rm = TRUE)
                app$normalizeGC$getBinZScore(sourceId, sample_name, binCounts, bd$bins$genome$pct_gc, bd$bins$genome$nAlleles)
            }, scoreType$distUnit)
            aggregateScores <- aggregateSampleScores(bd, sd$gcLimits, sd$stageTypes, sampleScores, scoreType$distUnit)
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
getScoreLevel <- function(scoreType){
         if(scoreType %in% names(scoreTypes$genome)) 'genome'
    else if(scoreType %in% names(scoreTypes$sample)) 'sample'
    else "NA"
}
getScoreType <- function(sourceId, scoreType){
    scoreLevel <- getScoreLevel(scoreType)
    paScores(sourceId)$scoreTypes[[scoreLevel]][[scoreType]]
}
getTypedStages <- function(sourceId) unlist(paScores(sourceId)$stageTypes)
getStageTypesByStage <- function(sourceId, stages) unlist(paScores(sourceId)$reverseStageTypes[stages])


#----------------------------------------------------------------------
# sample-level score structure summary and associated score object retrieval
#----------------------------------------------------------------------
# scores <- paScores(sourceId)$scores[[scoreLevel]]
# d print(names(scores))
# d print(names(scores[[input$scoreType]]))
# d print(names(scores[[input$scoreType]]$sampleScores)) # names are sample names
# d print(names(scores[[input$scoreType]]$aggregateScores))
# d print(names(scores[[input$scoreType]]$aggregateScores$by_stage)) # names are stage names
# d print(names(scores[[input$scoreType]]$aggregateScores$by_stageType)) # names are stage type names
# d print(names(scores[[input$scoreType]]$aggregateScores$stageType_delta)) # this is a scores object
# [1] "cpm"  "gcrz" "snif" "nrll"
# [1] "sampleScores"    "aggregateScores"
#  [1] "day35-wt-rs-rep1" "day35-wt-rs-rep2" "day38-wt-rs-rep1" "day38-wt-rs-rep2"
#  [5] "day32-wt-es"      "day34-wt-es"      "day35-wt-es-rep1" "day35-wt-es-rep2"
#  [9] "day38-wt-es-rep1" "day38-wt-es-rep2"
# [1] "by_stage"        "by_stageType"    "stageType_delta"
# [1] "early_round" "late_round"  "early_elong" "int_elong"   "late_elong"        
# [1] "round" "elong"
# [1] "score"    "dist"     "mean"     "sd"       "z"        "quantile"
getGenomeScores <- function(sourceId, scoreType){ # returns a single genome-level score object
    paScores(sourceId)$scores$genome[[scoreType]]
}
getSampleScoresList <- function(sourceId, scoreType){ # returns a list of sample-level score objects based on GC normalization
    if(scoreType == "gcrz"){
        gcResidualsZ(sourceId)
    } else {
        paScores(sourceId)$scores$sample[[scoreType]]
    }
}
getSampleScores <- function(sourceId, scoreType, samples){ # returns a list of sample-level score objects
    getSampleScoresList(sourceId, scoreType)$sampleScores[samples$sample_name]
}
getStageScores <- function(sourceId, scoreType, samples){ # returns a list of stage-level score objects matching a list of samples
    getSampleScoresList(sourceId, scoreType)$aggregateScores$by_stage[unique(samples$stage)]
}
getStageTypeScores <- function(sourceId, scoreType, samples){ # returns a list of stageType-level score objects matching a list of samples
    stageTypes <- getStageTypesByStage(sourceId, samples$stage)
    getSampleScoresList(sourceId, scoreType)$aggregateScores$by_stageType[unique(stageTypes)]
}
getStageTypeDeltaScores <- function(sourceId, scoreType){ # returns a single score object for the delta between stage types (or a genome-level score object)
    list(
        stageType_delta = if(getScoreLevel(scoreType) == "sample") {
            getSampleScoresList(sourceId, scoreType)$aggregateScores$stageType_delta
        } else {
            paScores(sourceId)$scores$genome[[scoreType]]
        }
    )
}
