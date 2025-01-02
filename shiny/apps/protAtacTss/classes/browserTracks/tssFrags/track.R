#----------------------------------------------------------------------
# tssFrags trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
# tssFragsBuffers <- reactiveValues()
# tfExpandReactive <- reactiveVal(NULL)
# tssFragsYBreaks <- data.table(height = numeric(), scoreTypeName = character(), rowType = character(), seriesName = character())

# constructor for the S3 class; REQUIRED
new_tssFragsTrack <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = TRUE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}

# build method for the S3 class; REQUIRED
build.tssFragsTrack <- function(track, reference, coord, layout){

    startSpinner(session, message = "getting tss data")
    sourceId <- track$settings$items()[[1]]$Source_ID
    req(sourceId)
    tfd <- paTssFragsData(sourceId)
    req(tfd)
    txnState <- track$settings$get("TSS_Frags","Transcription_State")
    tssFrags <- do.call(rbind, lapply(names(tfd$tssFrags[[txnState]]), function(sample_name_){
        tf <- tfd$tssFrags[[txnState]][[sample_name_]][chrom == coord$chromosome & start0 < coord$end & end1 >= coord$start]
        tf[, .(
            stage = tfd$samples[sample_name == sample_name_, stage],
            start1 = start0 + 1,
            end1   = end1,
            size   = end1 - start0
        )]
    }))
    req(nrow(tssFrags) > 0)

    stages <- tfd$samples[, unique(stage)]
    nStages <- length(stages)

    maxSize <- 600

    # set the plot frame
    padding <- padding(track, layout)
    height <- height(track, 2) + padding$total # or set a known, fixed height in inches
    ylim <- c(0, nStages * maxSize + 1)

    # make the plot
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # plot spans as lines
        color <- rgb(0, 0, 1, 0.25)
        I <- 1:nStages
        for(i in I){ 
            tf <- tssFrags[stage == stages[i]]  
            yOffset <- maxSize * (nStages - i)
            y <- tf$size + yOffset
            segments(tf$start1, y, tf$end1, y, col = color, lwd = 1.5)
            text(coord$start, yOffset + 50, stages[i], pos = 4, cex = 0.95)
            abline(h = yOffset, col = "black")
            abline(h = yOffset + 146, col = "grey")
        }
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

# track interaction methods
click.tssFragsTrack <- function(track, click, regionI){
    req(click$coord$y > 0)
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.tssFragsTrack <- showTrackSourcesDialog

# method for the S3 class to populate one or more trackNav inputs above the browser output
# only one navigation set is shown per track, your navigation should decide how to handle multiple regions
navigation.tssFragsTrack <- function(track, session, id, browser){

    # initialize the trackNavs, including input observers to handle user actions
    # initTrackNav will fail silenty if setting Track/Show_Navigation is not set or =="hide"
    navName <- initTrackNav(track, session, "activeTSS") # table reactive functions are provided below

    # as needed, create a reactive with the data to enumerate in the trackNavs
    trackNavData <- reactive({
        sourceId <- track$settings$items()[[1]]$Source_ID
        req(sourceId)
        txnState <- track$settings$get("TSS_Frags","Transcription_State")
        paTssFragsData(sourceId)$tss[[txnState]]
    })
    tagList(
        trackNavTable(
            track, 
            session, 
            browser$id,
            navName, # the name as provided by initTrackNav
            tableData = trackNavData, # populate a table based on track settings, etc.
            actionFn = function(selectedRow){
                req(selectedRow)
                d <- trackNavData()[selectedRow]
                handleTrackNavTableClick(1, track, d$chrom, d$start, d$end)
            }
        )
    )
}
