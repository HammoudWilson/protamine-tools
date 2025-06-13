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
getGenomeScores <- function(sourceId, scoreTypeName, binI){ # returns a single genome-level score object
    paScores_bins(sourceId)[binI][[scoreTypeName]]
}
getSampleScores <- function(metadata, config, coord, scoreTypeName, dataType){ # returns a list of sample-level score objects based on GC normalization
    scoresDir <- trimws(config$Scores_Dir)
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
    paScores_getCached(
        "bed_region_scores",
        keyObject = list(config$sourceId, config$bedFileName, scoreTypeName, dataType),
        from = 'disk',
        createFn = function(...) {
            scoresDir <- trimws(config$Scores_Dir)
            if(scoresDir == "") scoresDir <- metadata$env$SCORES_DIR
            scoreTable <- metadata$scoreTables[[scoreTypeName]][[dataType]]
            bgzFileName <- scoreTable$bgzFile
            bgzFile <- file.path(scoresDir, bgzFileName)
            req(file.exists(bgzFile))
            bins <- paScores_bins(config$sourceId) # remember, score bins are primary genome only
            I <- getIncludedAutosomeBins_scores(bins)
            d <- fread(
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
            )[I]
            overlaps <- foverlaps(
                d,                # larger table, smaller intervals
                config$bedData,   # keyed, smaller table, larger intervals
                type = "any",     # detect any overlap
                mult = "first",   # we only need to know if there is at least one overlap
                nomatch = NA,     # keep non-matching rows
                which = TRUE      # return indices instead of joined data
            )
            d[, hasOverlap := !is.na(overlaps)] 
            d
        },
        spinnerMessage = paste("parsing overlaps")
    )
}
getSeriesAggNames <- function(metadata, config){
    switch(
        config$Aggregate_By,
        sample     = metadata$samples$sample_name,
        stage      = unique(metadata$samples$stage),
        stage_type = names(metadata$stageTypes)
    )
}

#----------------------------------------------------------------------
# support score data download
#----------------------------------------------------------------------
getSampleScores_all <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects based on GC normalization
    metadata <- paScores_metadata(sourceId)
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
