#----------------------------------------------------------------------
# atacFootprint trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
# atacFootprintBuffers <- reactiveValues()
# atacFootprintExpandReactive <- reactiveVal(NULL)

# constructor for the S3 class; REQUIRED
new_atacFootprintTrack <- function(trackId) {
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
build.atacFootprintTrack <- function(track, reference, coord, layout){
    Max_Width_Bases <- track$settings$get("Footprint","Max_Width_Bases")
    if(!isTruthy(coord$width <= Max_Width_Bases)) return(trackInfo(track, coord, layout, "window too wide to plot footprint"))
    sourceId <- track$settings$items()[[1]]$Source_ID
    req(sourceId)
    metadata <- paTss_footprint(sourceId)
    req(metadata)
    
    # calculate plot parameters and dimensions to establish rules for binning and rendering
    plotWidthPixels <- as.integer(layout$plotWidth * layout$dpi)
    basesPerPixel <- coord$width / plotWidthPixels
    coordStart1 <- as.integer(coord$start) # comes in as bit64, incompatible with some functions below
    coordEnd1   <- as.integer(coord$end)
    pixelStart1s <- as.integer(coordStart1 + (1:plotWidthPixels - 1L) * basesPerPixel)
    pixelEnd1s <- c(pixelStart1s[-1] - 1L, coordEnd1)
    config <- list(
        sourceId        = sourceId,
        plotWidthPixels = plotWidthPixels,
        basesPerPixel   = basesPerPixel,
        coordStart1     = coordStart1,
        coordEnd1       = coordEnd1,
        Inserts_Bgz_Dir = track$settings$get("Data_Path","Inserts_Bgz_Dir")
    )
    config <- c(config, sapply(
        names(track$settings$Footprint()), 
        function(x) track$settings$get("Footprint", x),
        simplify = FALSE, USE.NAMES = TRUE
    ))
    startSpinner(session, message = "loading inserts")
    inserts <- paTss_get_inserts(metadata, coord, config)
    nSeries <- length(inserts)
    seriesNames <- names(inserts)

    # set the plot frame
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches
    ylim_series <- if(config$Plot_Inserts_As == "endpoint_counts"){
        # at present counts are normalized within the window, not to the genome
        y <- unlist(sapply(inserts, function(d) c(d$starts$y, d$ends$y)))
        c(0, 1.05 * quantile(y[y > 0], config$Endpoint_Max_Quantile))
    } else {
        c(0, config$Max_Insert_Size * 1.05)
    }
    ylim <- c(0, nSeries * 1)
    seriesRange <- (ylim_series[2] - ylim_series[1])

    # make the plot
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"
        for(i in 1:nSeries){ 
            yOffset <- (nSeries - i)
            abline(h = yOffset, col = "black")
            paTSS_add_series(inserts[[i]], yOffset, ylim_series, seriesRange, seriesNames[i], metadata, config)
            text(coord$start, yOffset + 0.8, seriesNames[i], pos = 4, cex = 1.25)
        }
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

# track interaction methods
click.atacFootprintTrack <- function(track, click, regionI){
    req(click$coord$y > 0)
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.atacFootprintTrack <- showTrackSourcesDialog
