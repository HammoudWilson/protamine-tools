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
        navigation = TRUE, # whether the track offers a custom, additional row of within-track navigation inputs
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
    config <- c(config, sapply(
        names(track$settings$Nucleosome_Chains()), 
        function(x) track$settings$get("Nucleosome_Chains", x),
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
    } else if(config$Plot_Inserts_As == "nucleosome_bias"){
        c(-1, 1)
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

        # overplot the called nucleosome chains as rectangles
        if(config$Show_Nucleosome_Chains && config$Aggregate_By == "stage"){
            regions <- paTss_ab_initio(sourceId)$regions[
                chrom  == coord$chrom & 
                start0 <  coordEnd1 & # wider than the plotted spans, includes the analysis flanks
                end1   >= coordStart1
            ]
            if(nrow(regions) > 0){

                # overplot the called overlap groups across all stages
                regions_group <- regions[index_stage == "overlap_group"]
                nGroups <- nrow(regions_group)
                if(nGroups > 0) rect(
                    regions_group$start0 + 1, 
                    rep(0, nGroups),
                    regions_group$end1, 
                    rep(nSeries, nGroups), 
                    border = CONSTANTS$plotlyColors$grey,
                    lwd = config$Overlay_Line_Width
                )

                # overplot the called nucleosome chains by stage
                for(stageI in 1:nSeries){
                    yOffset <- nSeries - stageI
                    regions_stage <- regions[index_stage == seriesNames[stageI]]
                    nChains <- nrow(regions_stage)
                    if(nChains == 0) next
                    rect(
                        regions_stage$start0 + 1, 
                        rep(yOffset + 0.025, nChains),
                        regions_stage$end1, 
                        rep(yOffset + 0.975, nChains), 
                        border = CONSTANTS$plotlyColors$red,
                        lwd = config$Overlay_Line_Width
                    )
                    for(chainI in 1:nChains){
                        nuc_starts0 <- as.integer(unlist(strsplit(regions_stage$nuc_starts0[chainI], ",")))
                        nNucs <- length(nuc_starts0)
                        rect(
                            nuc_starts0 + 1, 
                            rep(yOffset + 0.025, nNucs),
                            nuc_starts0 + 147, 
                            rep(yOffset + 0.975, nNucs),
                            border = NA, #CONSTANTS$plotlyColors$red,
                            col = CONSTANTS$plotlyColors$red %>% addAlphaToColor(0.1),
                            lwd = config$Overlay_Line_Width + 0.5
                        )
                    }
                }
            }
        }

        # plot the inserts
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

# method for the S3 class to populate one or more trackNav inputs above the browser output
# only one navigation set is shown per track, your navigation should decide how to handle multiple regions
navigation.atacFootprintTrack <- function(track, session, id, browser){
    Nav_Type <- track$settings$get("Navigation", "Nav_Type")

    # initialize the trackNavs, including input observers to handle user actions
    # initTrackNav will fail silenty if setting Track/Show_Navigation is not set or =="hide"
    navName <- initTrackNav(track, session, "dinuc_regions") # table reactive functions are provided below

    # as needed, create a reactive with the data to enumerate in the trackNavs
    if(Nav_Type == "plot"){
        getNavPlotRange <- function(d) {
            lim <- range(d, na.rm = TRUE)
            lim[1] <- lim[1] - diff(lim) * 0.05 # extend range by 5% on each side
            lim[2] <- lim[2] + diff(lim) * 0.05
            lim
        }
        trackNavPlotData <- reactive({
            sourceId <- track$settings$items()[[1]]$Source_ID
            req(sourceId)
            regions <- paTss_dinuc_regions(sourceId, "overlap_group")
            Plot_Nav_Type <- track$settings$get("Navigation", "Plot_Nav_Type")
            xy <- if(Plot_Nav_Type == "umap") {
                plotArgs <- list(
                    cex = 0.5,
                    xlab = "UMAP 1",
                    ylab = "UMAP 2"
                )
                as.matrix(regions[order(stage_mean), .(scaled_correlation_umap1, scaled_correlation_umap2)])
            } else {
                plotArgs <- list(
                    cex = 0.75,
                    xlab = "stage_mean",
                    ylab = "delta_RPKM"
                )
                as.matrix(regions[order(stage_mean), .(stage_mean, delta_RPKM)])
            }
            nRegions <- nrow(regions)
            color <- rainbow(nRegions / 0.9)[1:nRegions]
            xlim <- getNavPlotRange(xy[, 1])
            ylim <- getNavPlotRange(xy[, 2])
            I <- sample(1:nrow(xy))
            list(
                plotArgs = c(list(
                    x = xy[I, ],
                    col = color[I],
                    pch = 19
                ), plotArgs),
                layout = list(
                    width  = 96 * 6,
                    height = 96 * 3,
                    dpi = 96,
                    pointsize = 7,
                    mai = c(0.5, 0.5, 0.1, 0.1),
                    xlim = xlim,
                    ylim = ylim
                ),
                regions = regions[order(stage_mean)][I, .(chrom, start0, end1)]
            )
        })
        tagList(
            trackNavPlot(
                track, 
                session, 
                navName,
                contents = trackNavPlotData,
                actionFn = function(coord){
                    d <- trackNavPlotData()
                    req(coord, d)
                    dists <- rowSums((d$plotArgs$x - matrix(c(coord$x, coord$y), nrow(d$plotArgs$x), 2, byrow=TRUE))^2)
                    rowI <- which.min(dists)
                    region = d$regions[rowI]
                    handleTrackNavPlotClick(1, track, region$chrom, region$start0 + 1, region$end1)
                }
            )
        )
    } else {
        trackNavTableData <- reactive({
            sourceId <- track$settings$items()[[1]]$Source_ID
            req(sourceId)
            ai <- paTss_ab_initio(sourceId)
            Min_Max_RPKM <- track$settings$get("Navigation", "Table_Min_Max_RPKM") # filter the table a bit for faster loading
            Chromosome   <- track$settings$get("Navigation", "Table_Chromosome") %>% trimws()
            regions <- paTss_dinuc_regions(sourceId, "overlap_group")[max_RPKM >= Min_Max_RPKM]
            if (Chromosome != "") regions <- regions[chrom == Chromosome]
            paTss_parseRegionsForTable(sourceId, regions)
        })
        tagList(
            trackNavTable(
                track, 
                session, 
                browser$id,
                navName, # the name as provided by initTrackNav
                tableData = trackNavTableData, # populate a table based on track settings, etc.
                actionFn = function(selectedRow){
                    req(selectedRow)
                    d <- trackNavTableData()[selectedRow]
                    handleTrackNavTableClick(1, track, d$chrom, d$start0 + 1, d$end1)
                }
            )
        )
    }
}
