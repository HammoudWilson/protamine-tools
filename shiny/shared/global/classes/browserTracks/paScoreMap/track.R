#----------------------------------------------------------------------
# paScoreMap trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
paScoreBuffers    <- reactiveValues()
paExpandReactive  <- reactiveVal(NULL)
paScoreMapYBreaks <- data.table(height = numeric(), scoreTypeName = character(), rowType = character(), seriesName = character())

# constructor for the S3 class; REQUIRED
new_paScoreMapTrack <- function(trackId) {
    list(
        click = TRUE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = TRUE,
        items = TRUE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = TRUE,
        expand2 = FALSE
    )
}

# build method for the S3 class; REQUIRED
build.paScoreMapTrack <- function(track, reference, coord, layout){

    # calculate plot parameters and dimensions to establish rules for binning and rendering
    plotWidthPixels <- as.integer(layout$plotWidth * layout$dpi)
    basesPerPixel <- coord$width / plotWidthPixels
    coordStart1 <- as.integer(coord$start) # comes in as bit64, incompatible with some functions below
    coordEnd1   <- as.integer(coord$end)
    pixelStart1s <- as.integer(coordStart1 + (1:plotWidthPixels - 1L) * basesPerPixel)
    pixelEnd1s <- c(pixelStart1s[-1] - 1L, coordEnd1)
    config <- list(
        plotWidthPixels    = plotWidthPixels,
        sepHeightPixels    = 1L,
        headerHeightPixels = 16L,
        labelWidthPixels   = as.integer(layout$mai$left  * layout$dpi),
        legendWidthPixels  = as.integer(layout$mai$right * layout$dpi),
        Scores_Dir         = track$settings$get("Data_Path","Scores_Dir"),
        CutTag_Dir         = track$settings$get("Data_Path","CutTag_Dir"),
        basesPerPixel      = basesPerPixel,
        printMultiplier    = layout$printMultiplier
    )
    config <- c(config, sapply(
        names(track$settings$Score_Map()), 
        function(x) track$settings$get("Score_Map", x),
        simplify = FALSE, USE.NAMES = TRUE
    ))
    config$totalWidthPixels <- plotWidthPixels + config$labelWidthPixels + config$legendWidthPixels 

    # collect bins within the coordinate range
    startSpinner(session, message = "getting bins")
    sourceId <- track$settings$items()[[1]]$Source_ID
    req(sourceId)
    atac_metadata   <- paScores_metadata(sourceId)
    cuttag_metadata <- paCutTag_metadata(sourceId)
    bins <- paScores_bins(sourceId)
    req(bins)
    binI <- bins[chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start, binI]
    req(any(binI))
    config$samples <- paCutTag_samples(sourceId)
    config$binSize <- atac_metadata$env$BIN_SIZE
    config$binsPerPixel <- basesPerPixel / atac_metadata$env$BIN_SIZE

    # correlate bin to pixel crossing and overlap
    startSpinner(session, message = "parsing bin pixels")
    b <- bins[binI][, .( 
        binI, 
        start0,
        end1,
        included     = included, # at wider windows (the ones with more bins), most bins reside in single pixels
        pixel        = pmax(1,               floor((start0 + 1L - coordStart1) / basesPerPixel) + 1L),
        endPixel     = pmin(plotWidthPixels, floor((end1        - coordStart1) / basesPerPixel) + 1L)
    )]
    b <- rbind(
        b[pixel == endPixel, .( # pass single-pixel bins as is
            binI, 
            included,
            pixel,
            basesInPixel = atac_metadata$env$BIN_SIZE
        )],
        b[pixel != endPixel, {
            pixels <- pixel:endPixel
            .(
                included,
                pixel = pixels,
                basesInPixel = pmin(end1, pixelEnd1s[pixels]) - pmax(start0 + 1L, pixelStart1s[pixels]) + 1L 
            )
        }, by = .(binI)]
    )

    # ensure that excluded bins cannot contribute to track display images
    b[included == 0, basesInPixel := 0] 

    # use calls to scoreMapGroupImage to generate the single composite track image
    startSpinner(session, message = "rendering image")
    paScoreMapYBreaks <<- data.table(height = numeric(), scoreTypeName = character(), rowType = character(), seriesName = character())
    pngFile <- file.path(sessionDirectory, paste("paScoreMapTrack", "png", sep = "."))
    imager::imappend(
        lapply(
            # c("gc_z","txn","hic","stgm","gcrz_obs","nrll"), # ,"gcrz_wgt"
            config$Show_Tracks,
            scoreMapGroupImage,
            sourceId, atac_metadata, cuttag_metadata, config, coord, binI, b
        ),
        axis = 'y'
    ) %>% imager::save.image(pngFile)

    # save parameters for single-sample expansion 
    paScoreMapYBreaks[, y := cumsum(height)]
    maxY <- max(paScoreMapYBreaks$y)
    paScoreBuffers[[track$id]] <- list(
        coord = coord,
        sourceId = sourceId,
        config = config,
        b = b,
        maxY = maxY,
        binI = binI
    )

    # commit the final image
    list(
        ylim  = c(0, maxY),
        mai   = setMdiTrackMai(layout, padding(track, layout), mar = list(top = 0, bottom = 0)),
        image = pngToMdiTrackImage( # for tracks that generate images, not plots
            pngFile, 
            layout, 
            verticalPadding = 0L, # in pixels
            hasLeft  = TRUE,
            hasRight = TRUE
        )
    )
}

# track interaction methods
click.paScoreMapTrack <- function(track, click, regionI){
    req(click$coord$y > 0)
    app$browser$expandingTrack(regionI, list(trackId = track$id, row = row) )
}
brush.paScoreMapTrack <- function(track, brush, regionI){
    d <- paScoreBuffers[[track$id]]
    if(brush$keys$shift){
        showUserDialog(
            "Set as Training Region?", 
            tags$p("Would you like to set the selected region as a training region?"), 
            callback = function(...) {
                start0 <- brush$coord$x1 - 1
                end1   <- brush$coord$x2
                app$segmentation$addTrainingRegion(
                    d$sourceId, 
                    d$coord$chromosome,
                    start0, 
                    end1,
                    round((end1 - start0) / 1000, 1)
                )
            },
            fade = FALSE
        )
    } else {
        app$browser$jumpToCoordinates(
            regionI,
            d$coord$chromosome, 
            getX(brush$coord$x1), 
            getX(brush$coord$x2), 
            strict = TRUE
        )  
    }
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.paScoreMapTrack <- showTrackSourcesDialog

# expand method for the S3 class
# one expansion image can be shown per region, with same width as the main plots
# regionI must be passed to app$browser$expandingTrack
expand.paScoreMapTrack <- function(track, reference, coord, layout, regionI){
    config <- paScoreBuffers[[track$id]]$config
    req(config)
    startSpinner(session, message = "rendering color scales")
    padding <- padding(track, layout)

    n_scales <- 4

# z_score_color <- function(zScore, config, saturate = FALSE){
# quantile_score_color <- function(quantile, config, saturate = FALSE){
# txn_score_color <- function(log10cpm, config){
# cuttag_score_color <- function(rpkm, config){

# fraction_score_color <- function(fraction, config){

    height <- 0.5 * n_scales # or set a known, fixed height in inches
    # ylim <- range(y, na.rm = TRUE)
    ylim <- c(0, 1) # placeholder

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = ylim, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "row$scoreTypeName", # yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"


        points(0.5, 0.5, pch = 16, col = rgb(0, 0, 0, 0.5))





    })

    stopSpinner(session)

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}
