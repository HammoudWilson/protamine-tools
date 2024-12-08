# functions to load, parse, and filter bin data
protAtacLoadCreate <- "once" # for convenience when debugging load functions

# handle bin exclusions, including at gc extremes
isExcludedBin <- function(excluded, pct_gc){
    excluded == 1 | pct_gc < gcLimits[1] | pct_gc > gcLimits[2]
}
isIncludedAutosomeBin <- function(excluded, pct_gc, nAlleles){
    !isExcludedBin(excluded, pct_gc) & nAlleles == 2
}

# load and format genome bins and read counts
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

            # calculate the center of the genome bin GC distribution
            startSpinner(session, message = "loading bin counts..")
            bd$gc <- sapply(refTypes, function(refType){
                bd$bins[[refType]][
                    isIncludedAutosomeBin(excluded, pct_gc, nAlleles),
                    .(
                        mean = mean(pct_gc, na.rm = TRUE),
                        sd   = sd(  pct_gc, na.rm = TRUE)
                    )
                ]
            }, simplify = FALSE, USE.NAMES = TRUE)

            # calculate the center of the RPBA distribution
            # where RPBA is the Read Count Per Bin (of size bd$bin_size) Per Allele
            startSpinner(session, message = "loading bin counts...")
            bd$center <- sapply(refTypes, function(refType){
                I <- bd$bins[[refType]][, isIncludedAutosomeBin(excluded, pct_gc, nAlleles)] 
                sapply(bd$samples$sample_name, function(sample){
                    sapply(names(centerTypes_isZScore), function(centerType){
                        switch(
                            centerType,
                            rpba = {
                                bc <- rowSums(bd$binCounts[[refType]][I, , sample], na.rm = TRUE)
                                mean(bc, na.rm = TRUE) / 2
                            },
                            subnucleosomal_fraction = {
                                bc <- bd$binCounts[[refType]][I, , sample]
                                mean(bc[, "subnucleosomal"] / rowSums(bc, na.rm = TRUE), na.rm = TRUE)
                            },
                            0
                        )
                    }, simplify = FALSE, USE.NAMES = TRUE)
                }, simplify = FALSE, USE.NAMES = TRUE)   
            }, simplify = FALSE, USE.NAMES = TRUE)

            # create a share bin identifier for merge operations
            startSpinner(session, message = "loading bin counts....")
            for(refType in refTypes){
                bd$bins[[refType]]$binI <- 1:nrow(bd$bins[[refType]])
            }

            # establish the color palette for the unnormalized (raw) bin data
            bd$fold_change_col <- function(foldChange, maxFold){
                minFold <- 1 / maxFold 
                lfc <- log(pmax(minFold, pmin(maxFold, foldChange)), base = maxFold) # ranges from -1 to 1
                I <- floor(nTrackMapColorsPerSide * abs(lfc)) + 1L
                ifelse(lfc < 0, trackMapColors$low[I], trackMapColors$high[I])
            }
            bd$z_score_col <- function(zScore, maxZScore){
                minZScore <- -maxZScore 
                z <- pmax(minZScore, pmin(maxZScore, zScore))
                I <- floor(nTrackMapColorsPerSide * abs(z) / maxZScore) + 1L
                ifelse(z < 0, trackMapColors$low[I], trackMapColors$high[I])
            }
            bd
        }
    )$value
}

# load and format insert size data
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
            startSpinner(session, message = "loading inserts.")
            isd <- readRDS(getSourceFilePath(sourceId, "insertSizes"))
            # TODO: additional post-processing?
            isd
        }
    )$value
}
