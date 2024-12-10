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

# analyze and aggregate distributions of different bin scores
# all scores are expected to be one a comparable scale between samples, including
#   cpm  = read Counts Per Million reads
#   gcrz = GC Residual Z-Score (excess or deficit of reads relative to GC peers)
#   snif = Subnucleosomal Insert Fraction (fraction of bin insert sizes < 146bp)
#   nrll = subnucleosomal vs. nucleosomal Normalized Relative Log Likelihood
analyzeScoreDist <- function(bins, scores, distUnit = 0.1){
    scores_wrk <- scores[getIncudedAutosomeBins(bins)] # only use good autosomal bins to analyze score distributions
    dist <- data.table(x = as.integer(round(scores_wrk / distUnit)) * distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / sum(y, na.rm = TRUE)]
    mu <- mean(scores_wrk, na.rm = TRUE)
    sd <- sd(  scores_wrk, na.rm = TRUE)
    list(
        score = scores,
        dist  = dist,
        mean  = mu,
        sd    = sd,
        z     = (scores - mu) / sd
    )
}
analyzeSampleScores <- function(bd, scoreFn, distUnit = 0.1){
    sapply(bd$samples$sample_name, function(sample_name){
        analyzeScoreDist(bd$bins$genome, scoreFn(sample_name), distUnit)
    }, simplify = FALSE, USE.NAMES = TRUE)  
}
aggregateAndAnalyzeScores <- function(bins, analyzedScores, sample_names, distUnit = 0.1){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        analyzedScores[[sample_name]]$score
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(bins, apply(x, 1, mean, na.rm = TRUE), distUnit) # TODO: this is the slowest step, can we speed it up?
}
aggregateSampleScores <- function(bd, analyzedScores, distUnit = 0.1){
    # by_stage <- sapply(allStages, function(stage_){
    #     sample_names <- bd$samples[stage == stage_, sample_name]
    #     aggregateAndAnalyzeScores(bd$bins$genome, analyzedScores, sample_names, distUnit)
    # }, simplify = FALSE, USE.NAMES = TRUE)
    by_stageType <- sapply(names(stageTypes), function(stageType){
        dmsg(stageType)
        sample_names <- bd$samples[stage %in% stageTypes[[stageType]], sample_name]
        dstr(sample_names)
        aggregateAndAnalyzeScores(bd$bins$genome, analyzedScores, sample_names, distUnit)
    }, simplify = FALSE, USE.NAMES = TRUE)
    dstr(by_stageType)
    stageType_delta <- analyzeScoreDist(bd$bins$genome, by_stageType$elong$score - by_stageType$round$score, distUnit)
    dstr(stageType_delta)
    list(
        # by_stage = by_stage,
        by_stageType = by_stageType,
        stageType_delta = stageType_delta
    )
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

            # analyze the genome bin fractionGC distribution (also called pct_gc c/w bedtools nuc)
            startSpinner(session, message = "loading bin counts..")
            bd$gc <- sapply(refTypes, function(refType){
                analyzeScoreDist(bd$bins[[refType]], bd$bins[[refType]]$pct_gc, distUnit = 0.01)
                # bd$bins[[refType]][
                #     isIncludedAutosomeBin(excluded, pct_gc, nAlleles),
                #     .(
                #         mean = mean(pct_gc, na.rm = TRUE),
                #         sd   = sd(  pct_gc, na.rm = TRUE)
                #     )
                # ]
            }, simplify = FALSE, USE.NAMES = TRUE)

            # # calculate the center of the RPBA distribution
            # # where RPBA is the Read Count Per Bin (of size bd$bin_size) Per Allele
            # startSpinner(session, message = "loading bin counts...")
            # bd$center <- sapply(refTypes, function(refType){
            #     I <- bd$bins[[refType]][, isIncludedAutosomeBin(excluded, pct_gc, nAlleles)] 
            #     sapply(bd$samples$sample_name, function(sample){
            #         sapply(names(centerTypes_isZScore), function(centerType){
            #             switch(
            #                 centerType,
            #                 rpba = {
            #                     bc <- rowSums(bd$binCounts[[refType]][I, , sample], na.rm = TRUE)
            #                     mean(bc, na.rm = TRUE) / 2
            #                 },
            #                 subnucleosomal_fraction = {
            #                     bc <- bd$binCounts[[refType]][I, , sample]
            #                     mean(bc[, "subnucleosomal"] / rowSums(bc, na.rm = TRUE), na.rm = TRUE)
            #                 },
            #                 0
            #             )
            #         }, simplify = FALSE, USE.NAMES = TRUE)
            #     }, simplify = FALSE, USE.NAMES = TRUE)   
            # }, simplify = FALSE, USE.NAMES = TRUE)

            # create a share bin identifier for visualization merge operations
            startSpinner(session, message = "loading bin counts....")
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
            # TODO: additional post-processing?
            isd
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
    gcZScoreDeltaCache$get(
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
