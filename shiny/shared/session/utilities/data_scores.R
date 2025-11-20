#----------------------------------------------------------------------
# data recovery from atac/collate data packages
#----------------------------------------------------------------------
paScores_create <- "asNeeded"
paScores_getCached <- function(..., create = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(create)) create <- paScores_create
    protaminerCache$get(..., permanent = TRUE, create = create)$value
}

#----------------------------------------------------------------------
# load collate scores components into RAM
#----------------------------------------------------------------------
# load score summaries and distributions into RAM
paScores_metadata <- function(sourceId){
    paScores_getCached(
        "score_summary",
        key = sourceId,
        from = 'ram',
        createFn = function(...) {
            x <- readRDS(getSourceFilePath(sourceId, "scores"))
            x$bins <- NULL
            x$reverseStageTypes <- {
                reversed <- list()
                for (name in names(x$stageTypes)) {
                    for (value in x$stageTypes[[name]]) {
                        reversed[[value]] <- name
                    }
                }
                reversed
            }
            x
        },
        spinnerMessage = paste("loading score summary")
    )
}

# load score summaries and distributions into RAM
paScores_bins <- function(sourceId){
    paScores_getCached(
        "score_bins",
        key = sourceId,
        from = 'ram',
        createFn = function(...) {
            x <- readRDS(getSourceFilePath(sourceId, "scores"))
            x$bins[, ":="(
                binI = 1:.N,
                end1 = start0 + x$env$BIN_SIZE
            )]
            x$bins
        },
        spinnerMessage = paste("loading score bins")
    )
}

#----------------------------------------------------------------------
# score levels and stages
#----------------------------------------------------------------------
getScoreLevel <- function(scoreTypeName){
         if(scoreTypeName %in% names(scoreTypes$genome)) 'genome'
    else if(scoreTypeName %in% names(scoreTypes$sample)) 'sample'
    else "NA"
}
getScoreType <- function(scoreTypeName){
    scoreLevel <- getScoreLevel(scoreTypeName)
    scoreTypes[[scoreLevel]][[scoreTypeName]]
}
getTypedStages <- function(sourceId) unlist(paScores_metadata(sourceId)$stageTypes)
getStageTypesByStage <- function(sourceId, stages) unlist(paScores_metadata(sourceId)$reverseStageTypes[stages])

#----------------------------------------------------------------------
# score metadata retrieval (not the actual scores)
#----------------------------------------------------------------------
getSampleMetadataList <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects based on GC normalization
    paScores_metadata(sourceId)$scores$sample[[scoreTypeName]]
}
getSampleMetadata <- function(sourceId, scoreTypeName, samples){ # returns a list of sample-level score objects
    getSampleMetadataList(sourceId, scoreTypeName)$sampleScores[samples$sample_name]
}
getStageMetadata <- function(sourceId, scoreTypeName, samples){ # returns a list of stage-level score objects matching a list of samples
    getSampleMetadataList(sourceId, scoreTypeName)$aggregateScores$by_stage[unique(samples$stage)]
}
getStageTypeMetadata <- function(sourceId, scoreTypeName, samples){ # returns a list of stageType-level score objects matching a list of samples
    stageTypes <- getStageTypesByStage(sourceId, samples$stage)
    getSampleMetadataList(sourceId, scoreTypeName)$aggregateScores$by_stageType[unique(stageTypes)]
}
getStageTypeDeltaMetadata <- function(sourceId, scoreTypeName, cleanDist = FALSE){ # returns a single score object for the delta between stage types (or a genome-level score object)
    list(
        stageType_delta = if(getScoreLevel(scoreTypeName) == "sample") {
            getSampleMetadataList(sourceId, scoreTypeName)$aggregateScores$delta$stageType
        } else { # support plotting genome-level scores in the scoreSummary delta panel 
            x <- paScores_metadata(sourceId)$scores$genome[[scoreTypeName]]
            if(scoreTypeName == "txn" && cleanDist) x$dist <- x$dist[x %between% scoreTypes$genome$txn$valueLim]
            x
        }
    )
}

#----------------------------------------------------------------------
# score retrieval in a genome window
#----------------------------------------------------------------------
stgm_ecdf <- list()
hic_ecdf <- list()
getGenomeScores <- function(sourceId, scoreTypeName, binI, stgmQuantile = TRUE, hicQuantile = TRUE){ # returns a single genome-level score object
    bins <- paScores_bins(sourceId)
    scores <- bins[[scoreTypeName]]
    if(scoreTypeName == "stgm" && stgmQuantile) {
        if(is.null(stgm_ecdf[[sourceId]])) {
            scores_wrk <- scores[bins$included == 1 & !is.na(scores)]
            stgm_ecdf[[sourceId]] <<- ecdf(scores_wrk)
        }
        1 - stgm_ecdf[[sourceId]](scores[binI]) # thus, red = low stage mean (high in early stages), blue = high stage mean (high in late stages)
        # stgm_ecdf[[sourceId]](scores[binI])
    } else if(scoreTypeName == "hic" && hicQuantile) {
        if(is.null(hic_ecdf[[sourceId]])) {
            scores_wrk <- scores[bins$included == 1 & !is.na(scores)]
            hic_ecdf[[sourceId]] <<- ecdf(scores_wrk)
        }
        hic_ecdf[[sourceId]](scores[binI])
    } else {
        scores[binI]
    }
}
getSampleScores <- function(metadata, config, coord, scoreTypeName, dataType, isCutTagScore = FALSE){ # returns a list of sample-level score objects based on GC normalization
    scoresDir <- if(isCutTagScore) trimws(config$CutTag_Dir) else trimws(config$Scores_Dir)
    if(scoresDir == "") scoresDir <- metadata$env$SCORES_DIR
    scoreTable <- metadata$scoreTables[[scoreTypeName]][[dataType]]
    bgzFileName <- scoreTable$bgzFile
    bgzFile <- file.path(scoresDir, bgzFileName)
    req(file.exists(bgzFile))
    tabix <- getCachedTabix(bgzFile)
    getTabixRangeData(
        tabix, 
        coord, 
        col.names = scoreTable$colNames, 
        colClasses = c(
            "character", # chrom
            "integer",  # start0
            "integer",  # end1
            rep("numeric", length(scoreTable$colNames) - 3)
        )
    )
}
getSampleScores_regions <- function(metadata, config, scoreTypeName, dataType){ # returns a list of sample-level score objects based on GC normalization
    startSpinner(session, message = "getting bin overlaps")
    paScores_getCached(
        "bed_region_scores",
        keyObject = list(config$sourceId, config$bedData, scoreTypeName, dataType), # using bedData ensure uniqueness if file is recreated
        from = 'disk',
        createFn = function(...) {
            scoresDir <- trimws(config$Scores_Dir)
            if(scoresDir == "") scoresDir <- metadata$env$SCORES_DIR
            scoreTable <- metadata$scoreTables[[scoreTypeName]][[dataType]]
            bgzFileName <- scoreTable$bgzFile
            bgzFile <- file.path(scoresDir, bgzFileName)
            req(file.exists(bgzFile))
            bins <- paScores_bins(config$sourceId) # remember, score bins are primary genome only
            included_I <- getIncludedAutosomeBins_scores(bins)
            included_bin_scores <- fread(
                bgzFile, 
                header = FALSE, 
                sep = "\t",
                col.names = scoreTable$colNames, 
                colClasses = c(
                    "character", # chrom
                    "integer",   # start0
                    "integer",   # end1
                    rep("numeric", length(scoreTable$colNames) - 3)
                )
            )[included_I]
            overlaps <- foverlaps(
                x = included_bin_scores, # larger table, smaller intervals
                y = config$bedData,      # keyed, smaller table, larger intervals
                type = "any",            # detect any overlap
                mult = "first",          # we only need to know if there is at least one overlap
                nomatch = NA,            # keep non-matching rows from x, a.k.a, i, i.e., the first data.table
                which = TRUE             # return indices instead of joined data
            )
            # return all included bins with scores and a flag wether the bin overlapped any query region in bedData
            included_bin_scores[, hasOverlap := !is.na(overlaps)] 
            included_bin_scores
        },
        spinnerMessage = paste("parsing bin overlaps")
    )
}
getSeriesAggNames <- function(metadata, config, samplesFilter = TRUE){
    switch(
        config$Aggregate_By,
        sample     = metadata$samples[samplesFilter]$sample_name,
        stage      = unique(metadata$samples[samplesFilter]$stage),
        stage_type = names(metadata$stageTypes)
    )
}

#----------------------------------------------------------------------
# support score data download
#----------------------------------------------------------------------
getSampleScores_all <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects based on GC normalization
    cuttag <- scoreTypes$sample[[scoreTypeName]]$cuttag
    metadata <- if(!is.null(cuttag) && cuttag) paCutTag_metadata(sourceId)
                                          else paScores_metadata(sourceId)
    scoresDir <- metadata$env$SCORES_DIR
    scoreTable <- metadata$scoreTables[[scoreTypeName]]$score
    bgzFileName <- scoreTable$bgzFile
    bgzFile <- file.path(scoresDir, bgzFileName)
    req(file.exists(bgzFile))
    fread(
        cmd = paste("zcat", bgzFile),
        header = FALSE, 
        sep = "\t",
        col.names = scoreTable$colNames, 
        colClasses = c(
            "character", # chrom
            "integer",   # start0
            "integer",   # end1
            rep("numeric", length(scoreTable$colNames) - 3)
        )
    )
}
getSampleScores_allBins <- function(sourceId, scoreTypeName, stage){
    paScores_getCached(
        "getSampleScores_allBins",
        keyObject = list(sourceId, scoreTypeName, stage),
        from = 'disk',
        createFn = function(...) {
            getSampleScores_all(sourceId, scoreTypeName)[[stage]]
        },
        spinnerMessage = paste("loading stage scores")
    )
}
getStageTypeDelta_allBins <- function(sourceId, scoreTypeName){
    paScores_getCached(
        "getStageTypeDelta_allBins",
        keyObject = list(sourceId, scoreTypeName),
        from = 'ram',
        createFn = function(...) {
            getSampleScores_all(sourceId, scoreTypeName)$stageType
        },
        spinnerMessage = paste("loading score delta")
    )
}



    # env = env[c(
    #     'PRIMARY_GENOME',
    #     'SPIKE_IN_GENOME',
    #     'GENOME',
    #     'MAPPABILITY_SIZE_LEVELS',
    #     'MIN_INSERT_SIZE',
    #     'MAX_INSERT_SIZE',
    #     'BIN_SIZE',
    #     'HISTONE_STAGE',
    #     'PROTAMINE_STAGE',
    #     'SCORES_DIR'
    # )],
    # samples     = collate$samples,
    # references  = collate$references,
    # stageTypes  = stageTypes,
    # gcLimits    = gcLimits,
    # scores      = scores, # summary statistics only, see bed.bgz files for sample scores, bins for genome scores
    # scoreTables = scoreTables,
    # bins = collate$bins[isPrimaryGenome, .(
    #     chrom    = sub(paste0("-", env$PRIMARY_GENOME), "", chrom),
    #     start0   = start0,
    #     included = included,
    #     gc       = gc,
    #     gc_z     = gc_z,
    #     txn      = txn
    # )]
