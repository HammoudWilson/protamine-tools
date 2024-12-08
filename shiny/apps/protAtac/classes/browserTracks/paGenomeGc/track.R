#----------------------------------------------------------------------
# paGenomeGc trackBrowser track (i.e., a browserTrack)
# handles most facets of multi-sample ATAC-seq enrichment
#----------------------------------------------------------------------

# constructor for the S3 class; REQUIRED
new_paGenomeGcTrack <- function(trackId) {
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
build.paGenomeGcTrack <- function(track, reference, coord, layout){
    startSpinner(session, message = "getting bin data")
    sourceId <- track$settings$items()[[1]]$Source_ID
    bd <- paBinData(sourceId)
    windowBinI <- bd$bins$genome[, chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start]
    req(any(windowBinI))

    # calculate plot parameters
    Height_Pixels <- track$settings$get("Bin_Display","Height_Pixels")
    Max_Z_Score   <- track$settings$get("Bin_Display","Max_Z_Score")
    pixelWidth <- as.integer(layout$plotWidth * layout$dpi)
    basesPerPixel <- coord$width / pixelWidth
    coordStart1 <- as.integer(coord$start) # comes in as bit64, incompatible with some functions below
    coordEnd1   <- as.integer(coord$end)
    pixelStart1s <- as.integer(coordStart1 + (1:pixelWidth - 1L) * basesPerPixel)
    pixelEnd1s <- c(pixelStart1s[-1] - 1L, coordEnd1)
    startSpinner(session, message = "subsetting bin data")
    b <- bd$bins$genome[windowBinI][, .(
        binI, 
        start0,
        end1,
        excluded     = excluded,
        pct_gc       = pct_gc,
        pixel        = pmax(1,          floor((start0 + 1L - coordStart1) / basesPerPixel) + 1L), # at wider windows (the ones with more bins), most bins reside in single pixels
        endPixel     = pmin(pixelWidth, floor((end1        - coordStart1) / basesPerPixel) + 1L)
    )]
    b <- rbind(
        b[pixel == endPixel, .( # pass single-pixel bins as is
            binI, 
            excluded,
            pct_gc,
            pixel,
            basesInPixel = bd$bin_size
        )],
        b[pixel != endPixel, {
            pixels <- pixel:endPixel
            .(
                excluded,
                pct_gc,
                pixel = pixels,
                basesInPixel = pmin(end1, pixelEnd1s[pixels]) - pmax(start0 + 1L, pixelStart1s[pixels]) + 1L 
            )
        }, by = .(binI)]
    )
    b[excluded == 1, basesInPixel := 0]

    startSpinner(session, message = "aggregating bin data")
    gcz <- b[, 
        .(gcz = {
            z <- (pct_gc - bd$gc$genome$mean) / bd$gc$genome$sd
            weighted.mean(z, basesInPixel, na.rm = TRUE) 
        }),
        keyby = .(pixel)
    ][, 
        gcz
    ]
    startSpinner(session, message = "assigning colors")
    trackMatrix <-  matrix(
        rep(bd$z_score_col(
            gcz,
            Max_Z_Score
        ), Height_Pixels),
        ncol = 1
    )
    startSpinner(session, message = "rendering image")
    pngFile <- file.path(sessionDirectory, paste("paGenomeGcTrack", "png", sep = "."))
    matrixToCImg(trackMatrix, pixelWidth, Height_Pixels) %>% imager::save.image(pngFile)

    # commit the final image
    list(
        ylim  = c(0, 1),
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

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.paGenomeGcTrack <- showTrackSourcesDialog
