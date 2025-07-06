#----------------------------------------------------------------------
# ctcfMotifs trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
# ctcfMotifsBuffers <- reactiveValues()
# ctcfMotifsExpandReactive <- reactiveVal(NULL)

# constructor for the S3 class; REQUIRED
new_ctcfMotifsTrack <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}

# build method for the S3 class; REQUIRED
build.ctcfMotifsTrack <- function(track, reference, coord, layout){
    Max_CTCF_Window_Bp <- track$settings$get("CTCF_Motifs","Max_CTCF_Window_Bp")
    if(!isTruthy(coord$width <= Max_CTCF_Window_Bp)) return(trackInfo(track, coord, layout, "window too wide to plot CTCF motifs"))
    ctcf <- paTss_ctcf_motifs(reference)[
        chrom == coord$chrom & 
        pos1  >= coord$start & 
        pos1  <= coord$end
    ]

    # set the plot frame
    padding <- padding(track, layout)
    height <- 0.25 + padding$total # or set a known, fixed height in inches
    ylim <- c(0, 2)

    # make the plot
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"
        points(
            x = ctcf$pos1, 
            y = rep(1, nrow(ctcf)),
            pch = ifelse(ctcf$strand == "+", ">", "<"), # triangle up for + strand, triangle down for - strand
            cex = 2, 
            col = CONSTANTS$plotlyColors$brown
        )
        mtext("CTCF", side = 2, line = 0.5, at = 1, las = 1)
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}
