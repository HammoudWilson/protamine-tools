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

            # TODO: remove this and constants.R object once pipeline is re-run
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
gcResidualZScores <- function(sourceId, gcBiasModels = NULL){
    if(is.null(gcBiasModels)) app$normalizeGC$getGcBiasModels_externalCall(sourceId)
    gcrzCache$get(
        "gcZScoreDelta",
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
            analyzedScores <- analyzeSampleScores(bd, function(sample_name){
                startSpinner(session, message = paste("loading gcrz", sample_name))
                binCounts <- rowSums(bd$binCounts$genome[, , sample_name], na.rm = TRUE)
                x <- app$normalizeGC$getBinZScore(sourceId, sample_name, binCounts, bd$bins$genome$pct_gc, bd$bins$genome$nAlleles)
            })
            startSpinner(session, message = "aggregating gcrz")
            x <- aggregateSampleScores(bd, analyzedScores)

            dstr(x)

            x
        }
    )$value
}

#----------------------------------------------------------------------
# sample-level score structure summary and associated score object retrieval
#----------------------------------------------------------------------
# scores <- paScores(sourceId)$scores[[scoreLevel]]
# dprint(names(scores))
# dprint(names(scores[[input$scoreType]]))
# dprint(names(scores[[input$scoreType]]$sampleScores)) # names are sample names
# dprint(names(scores[[input$scoreType]]$aggregateScores))
# dprint(names(scores[[input$scoreType]]$aggregateScores$by_stage)) # names are stage names
# dprint(names(scores[[input$scoreType]]$aggregateScores$by_stageType)) # names are stage type names
# dprint(names(scores[[input$scoreType]]$aggregateScores$stageType_delta)) # this is a scores object
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
getSampleScores <- function(sourceId, scoreType, samples){ # returns a list of sample-level score objects
    paScores(sourceId)$scores$sample[[scoreType]]$sampleScores[samples$sample_name]
}
getStageScores <- function(sourceId, scoreType, samples){ # returns a list of stage-level score objects matching a list of samples
    paScores(sourceId)$scores$sample[[scoreType]]$aggregateScores$by_stage[unique(samples$stage)]
}
getStageTypeScores <- function(sourceId, scoreType, samples){ # returns a list of stageType-level score objects matching a list of samples
    stageTypes <- getStageTypesByStage(sourceId, samples$stage)
    paScores(sourceId)$scores$sample[[scoreType]]$aggregateScores$by_stageType[unique(stageTypes)]
}
getStageTypeDeltaScores <- function(sourceId, scoreType){ # returns a single score object for the delta between stage types
    list(
        stageType_delta = if(getScoreLevel(scoreType) == "sample") {
            paScores(sourceId)$scores$sample[[scoreType]]$aggregateScores$stageType_delta
        } else {
            paScores(sourceId)$scores$genome[[scoreType]]
        }
    )
}
