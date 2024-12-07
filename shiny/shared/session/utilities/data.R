# functions to load, parse, and filter bin data
protAtacLoadCreate <- "once" # for convenience when debugging load functions

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
        permanent = FALSE, # don't resave to disk, is large and slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){

            # load the bin data from the pipeline data package
            startSpinner(session, message = "loading bin counts.")
            bd <- readRDS(getSourceFilePath(sourceId, "binCounts"))

            # calculate the center of the bin GC distribution
            startSpinner(session, message = "loading bin counts..")
            bd$gc_center <- sapply(refTypes, function(refType){
                bd$bins[[refType]][
                    excluded == 0,
                    mean(pct_gc, na.rm = TRUE)
                ]
            }, simplify = FALSE, USE.NAMES = TRUE)

            # calculate the center of the RPBA distribution
            # where RPBA is the Read Count Per Bin (of size bd$bin_size) per allele
            startSpinner(session, message = "loading bin counts...")
            bd$center <- sapply(refTypes, function(refType){
                I <- bd$bins[[refType]][, excluded == 0 & nAlleles == 2] # don't use sex chroms for centering
                sapply(bd$samples$sample_name, function(sample){
                    sapply(centerTypes, function(centerType){
                        switch(
                            centerType,
                            rpba = {
                                bc <- rowSums(bd$binCounts[[refType]][I, , sample], na.rm = TRUE)
                                mean(bc, na.rm = TRUE) / 2
                            },
                            subnucleosomal_fraction = {
                                bc <- bd$binCounts[[refType]][I, , sample]
                                mean(bc[, "subnucleosomal"] / rowSums(bc, na.rm = TRUE), na.rm = TRUE)
                            }
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
            bd$raw <- {
                colorsPerSide <- 30
                list(
                    pltLow = colorRampPalette(c(paColors$GREY, paColors$BLUE))(colorsPerSide + 1), # blue color is cold/depleted
                    pltHgh = colorRampPalette(c(paColors$GREY, paColors$RED))( colorsPerSide + 1), # red  color is hot/ enriched
                    col = function(bd, foldChange, maxFold){
                        minFold <- 1 / maxFold 

                        # calculate the log fold change, with bounds, in base maxFold
                        # thus, resulting values range from -1 to 1
                        lfc <- log(pmax(minFold, pmin(maxFold, foldChange)), base = maxFold)

                        # calculate the color index
                        # floor component ranges from 0 to colorsPerSide
                        # final value ranges from 1 to colorsPerSide + 1
                        I <- floor(colorsPerSide * abs(lfc)) + 1L

                        # return the final colors depending on which side of neutral the value lies
                        ifelse(lfc < 0, bd$raw$pltLow[I], bd$raw$pltHgh[I])
                    }
                )
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
        permanent = FALSE, # don't resave to disk, is large and slow to load that way
        from = "ram",
        create = protAtacLoadCreate,
        createFn = function(...){

            # load the bin data from the pipeline data package
            startSpinner(session, message = "loading inserts.")
            isd <- readRDS(getSourceFilePath(sourceId, "insertSizes"))
            isd
        }
    )$value
}
