#----------------------------------------------------------------------
# paSampleComposite trackBrowser track (i.e., a browserTrack)
# handles most facets of multi-sample ATAC-seq enrichment
#----------------------------------------------------------------------

# constructor for the S3 class; REQUIRED
new_paSampleCompositeTrack <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}

# build method for the S3 class; REQUIRED
build.paSampleCompositeTrack <- function(track, reference, coord, layout){

    startSpinner(session, message = "getting bin data")
    sourceId <- track$settings$items()[[1]]$Source_ID
    bd <- paBinData(sourceId)
    windowBinI <- bd$bins$genome[, chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start]
    req(any(windowBinI))

    # calculate plot parameters
    Sample_Height_Pixels <- track$settings$get("Bin_Display","Sample_Height_Pixels")
    Max_Axis_Value       <- track$settings$get("Bin_Display","Max_Axis_Value")
    Bin_Data_Type        <- track$settings$get("Bin_Display","Bin_Data_Type")
    sampleNames <- bd$samples$sample_name
    nSamples <- length(sampleNames)
    pixelWidth <- as.integer(layout$plotWidth * layout$dpi)
    basesPerPixel <- coord$width / pixelWidth
    coordStart1 <- as.integer(coord$start) # comes in as bit64, incompatible with some functions below
    coordEnd1   <- as.integer(coord$end)
    pixelStart1s <- as.integer(coordStart1 + (1:pixelWidth - 1L) * basesPerPixel)
    pixelEnd1s <- c(pixelStart1s[-1] - 1L, coordEnd1)

    startSpinner(session, message = "subsetting bin data")
    binI <- bd$bins$genome[windowBinI, binI]
    b <- bd$bins$genome[windowBinI][, .(
        binI, 
        start0,
        end1,
        excluded     = excluded,
        # pct_gc       = pct_gc,
        pixel        = pmax(1,          floor((start0 + 1L - coordStart1) / basesPerPixel) + 1L), # at wider windows (the ones with more bins), most bins reside in single pixels
        endPixel     = pmin(pixelWidth, floor((end1        - coordStart1) / basesPerPixel) + 1L)
    )]
    b <- rbind(
        b[pixel == endPixel, .( # pass single-pixel bins as is
            binI, 
            excluded,
            # pct_gc,
            pixel,
            basesInPixel = bd$bin_size
        )],
        b[pixel != endPixel, {
            pixels <- pixel:endPixel
            .(
                excluded,
                # pct_gc,
                pixel = pixels,
                basesInPixel = pmin(end1, pixelEnd1s[pixels]) - pmax(start0 + 1L, pixelStart1s[pixels]) + 1L 
            )
        }, by = .(binI)]
    )
    b[excluded == 1, basesInPixel := 0]

    startSpinner(session, message = "aggregating bin data")
    nAlleles <- bd$bins$genome[windowBinI, nAlleles[1]]
    d <- sapply(sampleNames, function(sampleName){
        x <- merge(
            b, 
            data.table(
                binI = binI, 
                val = {
                    bc <- bd$binCounts$genome[windowBinI, , sampleName]
                    switch(
                        Bin_Data_Type,
                        rpba = rowSums(bc),
                        gc_adj_z_score = {
                            gc <- bd$bins$genome[windowBinI, pct_gc]
                            app$normalizeGC$getBinZScore(sourceId, sampleName, rowSums(bc), gc, nAlleles)
                        },
                        subnucleosomal_fraction = bc[, "subnucleosomal"] / rowSums(bc, na.rm = TRUE), 
                        subnucleosomal_nrll = paSegmentation(sourceId)$NRLL[windowBinI][[sampleName]]
                    )
                }
            ),
            by = "binI", 
            all.x = TRUE
        )
        switch(
            Bin_Data_Type,
            rpba = x[, 
                .(val = {
                    sumCount <- sum(val * basesInPixel / bd$bin_size, na.rm = TRUE)
                    nBins    <- sum(      basesInPixel / bd$bin_size, na.rm = TRUE)
                    if(nBins == 0) NA_real_ else sumCount / nBins / nAlleles
                }),
                keyby = .(pixel)
            ],
            x[, 
                .(val = {
                    weighted.mean(val, basesInPixel, na.rm = TRUE)
                }),
                keyby = .(pixel)
            ]
        )$val
    }, simplify = FALSE, USE.NAMES = TRUE)

    startSpinner(session, message = "assigning colors")
    trackMatrix <-  sapply(sampleNames, function(sampleName){
        cols <- if(centerTypes_isZScore[[Bin_Data_Type]]) bd$z_score_col(
            d[[sampleName]],
            Max_Axis_Value
        ) else bd$fold_change_col(
            d[[sampleName]] / bd$center$genome[[sampleName]][[Bin_Data_Type]],
            Max_Axis_Value
        )
        rep(cols, Sample_Height_Pixels)
    })

    startSpinner(session, message = "rendering image")
    pngFile <- file.path(sessionDirectory, paste("paSampleCompositeTrack", "png", sep = "."))
    # nImgRows <- nValTypes * nSamples + nValTypes - 1
    trackImgHeight <- Sample_Height_Pixels * nSamples
    matrixToCImg(trackMatrix, pixelWidth, trackImgHeight) %>% imager::save.image(pngFile)
    # legend <- cnvColorLegend(track, nImgRows)
    # imager::imappend(
    #     list(
    #         matrixToCImg(trackMatrix,   pixelWidth,   trackImgHeight), 
    #         matrixToCImg(legend$matrix, legend$width, trackImgHeight)
    #     ),
    #     axis = 'x'
    # ) %>% imager::save.image(pngFile)

    # # save parameters for single-sample expansion
    # msvCompositeBuffers[[track$id]] <- list(
    #     coord = coord,
    #     sourceId = sourceId,
    #     sampleIds = rev(c(sampleIds, NA, sampleIds, NA, sampleIds, NA, sampleIds)),
    #     valueTypes = rev(c(rep("LRR", nSamples), NA, rep("ZYG", nSamples), NA, rep("high", nSamples), NA, rep("low", nSamples)))
    # )

    # commit the final image
    list(
        ylim  = c(0, nSamples),
        mai   = setMdiTrackMai(layout, padding(track, layout), mar = list(top = 0, bottom = 0)),
        image = pngToMdiTrackImage( # for tracks that generate images, not plots
            pngFile, 
            layout, 
            verticalPadding = 5L, # in pixels
            hasLeft  = FALSE,
            hasRight = FALSE
        )
    )
}

# # plot interaction methods for the S3 class
# # called by trackBrowser if track$click, $hover, or $brush is TRUE, above
# # regionI indicates the region plot the user interacted with, the value must be passed to app$browser$jumpToCoordinates, etc.
# setMsvZoomSample <- function(coord, sourceId, sampleId, trackId, resolution){
#     msvZoomSample(list( # the sample currently in view in msvSampleZoom
#         coord       = coord,
#         sourceId    = sourceId,
#         sampleId    = sampleId,
#         trackId     = trackId, 
#         resolution  = resolution
#     ))
#     app$browser$forceTrackTypeRefresh("msvSampleZoom")
# }
# setMsvWorkingSample <- function(trackId, sourceId, sampleId, chrom, start = NULL, end = NULL){
#     msvWorkingSample(data.table( # the sample for which a CNV is currently being constructed
#         sourceId    = sourceId,
#         sampleId    = sampleId,
#         trackId     = trackId
#     ))
#     msvWorkingBoundaries(data.table(
#         sampleId    = sampleId,
#         chrom       = chrom,
#         position    = if(is.null(start)) NA_integer_ else c(start, end)
#     ))
# }
# click.paSampleCompositeTrack <- function(track, click, regionI){
#     req(click$coord$y > 0)
#     d <- msvCompositeBuffers[[track$id]]
#     req(d$sampleIds)
#     sampleId  <- d$sampleIds[ceiling(click$coord$y)]
#     valueType <- d$valueTypes[ceiling(click$coord$y)]
#     req(sampleId, valueType)

#     # all composite track clicks set the sample currently in view in msvSampleZoom
#     setMsvZoomSample(
#         coord      = d$coord,
#         sourceId   = d$sourceId,
#         sampleId   = sampleId,
#         trackId    = track$id,
#         resolution = if(valueType == "low") "low" else "high"
#     )

#     # shift-click jumps to the selected CNV, if any
#     # ctrl-click sets the working sample and jumps to the selected CNV, if any
#     if(click$keys$shift || click$keys$ctrl){
#         hmmCnvs <- msv_getTrackCnvs(track, coord = d$coord, sampleId_ = sampleId, filterRes = FALSE)[resolution == valueType]
#         hmmCnv <- hmmCnvs[start <= click$coord$x & end >= click$coord$x]
#         if(nrow(hmmCnv) > 0){
#             if(click$keys$ctrl) setMsvWorkingSample (track$id, d$sourceId, sampleId, d$coord$chromosome, hmmCnv$start, hmmCnv$end)
#             app$browser$jumpToCoordinates(1, d$coord$chromosome, hmmCnv$start, hmmCnv$end)
#         } else if(click$keys$ctrl){
#             setMsvWorkingSample(track$id, d$sourceId, sampleId, d$coord$chromosome)
#         }
#     }
# }

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.paSampleCompositeTrack <- showTrackSourcesDialog

# # method for the S3 class to populate one or more trackNav inputs above the browser output
# # only one navigation set is shown per track, your navigation should decide how to handle multiple regions
# navigation.paSampleCompositeTrack <- function(track, session, id, browser){

#     # initialize the trackNavs, including input observers to handle user actions
#     msvCnvNavName <- initTrackNav(track, session, "msvCnvNavTable") # table reactive functions are provided below

#     # as needed, create a reactive with the data to enumerate in the trackNavs
#     msvCnvNavTable <- reactive({
#         inactivateParsedCnvs[[track$settings$items()[[1]]$Source_ID]]
#         msv_getTrackCnvs(track, filterRes = TRUE)[, .(
#             sampleId,
#             chrom,
#             start,
#             end,
#             size_kb,
#             cnvType = cnvType_,
#             CNC = CNC,
#             nProbes,
#             nInf,
#             RLL = ifelse(RLL == Inf, track$settings$get("MSV_Display","MaxDisplay_RLL"), RLL),
#             LRR,
#             ZYG, 
#             status
#         )]
#     })

#     # handle table click
#     clickHandler <- function(selectedRow){
#         req(selectedRow)
#         hmmCnv <- msvCnvNavTable()[selectedRow]
#         handleTrackNavTableClick(NULL, track, hmmCnv$chrom, hmmCnv$start, hmmCnv$end, expandFn)
#     }

#     # return the associated navigation UI elements
#     tagList(
#         trackNavTable(
#             track, 
#             session, 
#             browser$id,
#             msvCnvNavName, # the name as provided by initTrackNav
#             tableData = msvCnvNavTable, # populate a table based on track settings, etc.
#             actionFn = clickHandler # add other argument to pass to bufferedTableServer, but omit "selection" and "width"
#         )
#     )
# }
