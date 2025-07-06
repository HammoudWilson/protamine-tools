#----------------------------------------------------------------------
# server components for the dinucRegions appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
dinucRegionsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'dinucRegions'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
regions <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTss_dinuc_regions(sourceId)
})

#----------------------------------------------------------------------
# interactive umap plot
#----------------------------------------------------------------------
umapPlotData <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    ai <- paTss_ab_initio(sourceId)
    regions <- regions()
    req(ai, regions)
    startSpinner(session, message = "loading umap data")
    Dinuc_Scaling <- umapPlotBox$settings$get("UMAP", "Dinuc_Scaling")
    Dinuc_Metric  <- umapPlotBox$settings$get("UMAP", "Dinuc_Metric")
    nRegions <- nrow(regions)
    list(
        xy      = ai$umap[[Dinuc_Scaling]][[Dinuc_Metric]][order(regions$stage_mean), ],
        regions = regions[order(stage_mean)],
        title   = paste(Dinuc_Scaling, Dinuc_Metric, sep = " - "),
        color   = rainbow(nRegions / 0.9)[1:nRegions]
    )
})
createUmapPlot <- function(settings, plot){
    d <- umapPlotData()
    req(d)
    startSpinner(session, message = "rendering umap plot")
    par(mar = c(4, 4, 0, 0) + 0.1) 
    layout <- plot$initializePng() %>% plot$initializeFrame(
        xlim = range(d$xy[, 1]),
        ylim = range(d$xy[, 2]),
        xlab = "UMAP 1",
        ylab = "UMAP 2"
        ,
        # # cex.main = 0.95,
        # title = d$title
    )
    I <- sample(1:nrow(d$xy))
    plot$addPoints(
        x = d$xy[I, 1],
        y = d$xy[I, 2],
        col = d$color[I],
        pch = settings$get("Points_and_Lines", "Point_Type"),
        cex = settings$get("Points_and_Lines", "Point_Size")
    )
    stopSpinner(session)
    plot$finishPng(layout)
}
umapPlotBox <- mdiInteractivePlotBoxServer(
    "umapPlotBox",
    click = TRUE,
    points = TRUE,
    settings = list(
        UMAP = list(
            Dinuc_Scaling = list(
                type = "selectInput",
                choices = c("scaled", "unscaled"),
                value = "scaled"
            ),
            Dinuc_Metric  = list(
                type = "selectInput",
                choices = c("euclidean", "cosine", "correlation"),
                value = "correlation"
            )
        )
    ),
    defaults = list(
        Plot_Frame = list(
            Width_Inches = 4,
            Height_Inches = 3
        ),
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.25
        )
    ),
    create = function(...) createUmapPlot(..., umapPlotBox) # a function or reactive that creates the plot as a png file using settings and helpers
)
observeEvent(umapPlotBox$plot$click(), {
    setSelectedRegion(umapPlotBox$plot$click()$coord, umapPlotData())
})

#----------------------------------------------------------------------
# interactive correlation plot
#----------------------------------------------------------------------
correlationPlotData <- reactive({
    regions <- regions()
    req(regions)
    startSpinner(session, message = "loading correlation")
    regions <- regions[order(stage_mean)]
    nRegions <- nrow(regions)
    x_axis <- correlationPlotBox$settings$get("Correlation", "X_Axis")
    y_axis <- correlationPlotBox$settings$get("Correlation", "Y_Axis")
    list(
        x_axis = x_axis,
        y_axis = y_axis,
        xy = as.matrix(regions[, .SD, .SDcols = c(x_axis, y_axis)]),
        regions = regions,
        color = rainbow(nRegions / 0.9)[1:nRegions]
    )
})
getCorrelationRange <- function(d, col) {
    lim <- range(d, na.rm = TRUE)
    lim[1] <- lim[1] - diff(lim) * 0.05 # extend range by 5% on each side
    lim[2] <- lim[2] + diff(lim) * 0.05
    lim
}
createCorrelationPlot <- function(settings, plot){
    d <- correlationPlotData()
    req(d)
    startSpinner(session, message = "rendering correlation")
    par(mar = c(4, 4, 0, 0) + 0.1) 
    layout <- plot$initializePng() %>% plot$initializeFrame(
        xlim = getCorrelationRange(d$xy[, 1]),
        ylim = getCorrelationRange(d$xy[, 2]),
        xlab = d$x_axis,
        ylab = d$y_axis,
        xaxs = "i",
        yaxs = "i"
    )
    I <- sample(1:nrow(d$xy))
    plot$addPoints(
        x = d$xy[I, 1],
        y = d$xy[I, 2],
        col = d$color[I],
        pch = settings$get("Points_and_Lines", "Point_Type"),
        cex = settings$get("Points_and_Lines", "Point_Size")
    )
    stopSpinner(session)
    plot$finishPng(layout)
}
correlationPlotBox <- mdiInteractivePlotBoxServer(
    "correlationPlotBox",
    click = TRUE,
    points = TRUE,
    settings = list(
        Correlation = list(
            X_Axis = list(
                type = "selectInput",
                choices = paTss_appRPKMCols,
                value = "stage_mean"
            ),
            Y_Axis = list(
                type = "selectInput",
                choices = paTss_appRPKMCols,
                value = "delta_RPKM"
            )
        )
    ),
    defaults = list(
        Plot_Frame = list(
            Width_Inches = 4,
            Height_Inches = 3
        ),
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    ),
    create = function(...) createCorrelationPlot(..., correlationPlotBox)
)
observeEvent(correlationPlotBox$plot$click(), {
    setSelectedRegion(correlationPlotBox$plot$click()$coord, correlationPlotData())
})

#----------------------------------------------------------------------
# selected region profile plot
#----------------------------------------------------------------------
selectedRegion <- reactiveVal(NULL)
setSelectedRegion <- function(coord, d) {
    dists <- rowSums((d$xy - matrix(c(coord$x, coord$y), nrow(d$xy), 2, byrow=TRUE))^2)
    rowI <- which.min(dists)
    selectedRegion(list(
        region = d$regions[rowI],
        color  = d$color[rowI]
    ))
}
regionProfilePlot <- staticPlotBoxServer(
    "regionProfilePlot",
    title = TRUE,
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        ai <- paTss_ab_initio(sourceId)
        d <- selectedRegion()
        req(ai, d)
        startSpinner(session, message = "rendering profile")
        nStages <- length(ai$stages)
        regionProfilePlot$initializeFrame(
            xlim = c(1, nStages),
            ylim = c(0, d$region$max_RPKM * 1.05),
            xlab = "Stage",
            ylab = "RPKM",
            title = paste0(d$region$chrom, ":", d$region$start0, "-", d$region$end1)
        )
        regionProfilePlot$addLines(
            x = 1:nStages,
            y = unlist(d$region[, .SD, .SDcols = ai$stages]),
            col = d$color,
            lwd = 2
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# selected region insert footprint plot
#----------------------------------------------------------------------
createRegionPlot <- function(settings, plot) {
    region <- selectedRegion()$region
    req(region)
    startSpinner(session, message = "rendering region footprint")
    padding_bp <- 250
    coord <- list(chromosome = region$chrom, start = region$start0 + 1 - padding_bp, end = region$end1 + padding_bp)
    coord$range <- c(coord$start, coord$end)
    coord$width <- coord$end - coord$start + 1
    app$browser$createBrowserPlot(
        regionI = 1, 
        pngFile = plot$pngFile,
        externalCoord = coord
    )$layout
}
regionPlotBox <- mdiInteractivePlotBoxServer(
    "regionPlotBox",
    defaults = list(
        Plot_Frame = list(
            Width_Inches = 8,
            Height_Inches = 4
        )
    ),
    create = function(...) createRegionPlot(..., regionPlotBox) # a function or reactive that creates the plot as a png file using settings and helpers
)

#----------------------------------------------------------------------
# profile plot by cluster/quantile
#----------------------------------------------------------------------
clusterSettings <- list(
    Clustering = list(
        Cluster_Col = list(
            type = "selectInput",
            choices = c("quantile", "k_scaled", "k_unscaled"),
            value = "quantile"
        ),
        Plots_Clusters_As = list(
            type = "selectInput",
            choices = c("traces", "heatmap"),
            value = "traces"
        )
    )
)
clusterAggregates <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    ai <- paTss_ab_initio(sourceId)
    regions <- paTss_dinuc_regions(sourceId)
    req(ai, regions)
    startSpinner(session, message = "aggregating clusters")
    Cluster_Col <- clusterProfilePlot$settings$get("Clustering", "Cluster_Col")
    clusters <- sapply(ai$stages, function(stage) regions[, mean(regions[[stage]][.I]), keyby = Cluster_Col][[2]])
    clusterCounts <- regions[, .N, keyby = Cluster_Col]
    nClusters <- nrow(clusters)
    nStages   <- ncol(clusters)
    stage_means <- sapply(1:nClusters, function(i) weighted.mean(1:nStages, clusters[i,]))
    if(Cluster_Col != "quantile") {
        I <- order(stage_means)
        clusters <- clusters[I, ]
        stage_means <- stage_means[I]
        clusterCounts <- clusterCounts[I]
    }
    list(
        ai = ai,
        clusters = clusters,
        clusterCounts = clusterCounts,
        stage_means = stage_means,
        nClusters = nClusters,
        nStages = nStages,
        stages = ai$stages,
        colors = rainbow(nClusters / 0.9)[1:nClusters]
    )
})
plotClusterTraces <- function(d) {
    par(mar = c(4, 4, 0, 0) + 0.1)
    ylim <- c(0, max(d$clusters, na.rm = TRUE) * 1.05)
    clusterProfilePlot$initializeFrame(
        xlim = c(1, d$nStages),
        ylim = ylim,
        xlab = "",
        ylab = "Average RPKM",
        xaxt = "n"
    )
    axis(side = 1, labels = FALSE)
    text(1:d$nStages, -ylim[2] / 8, d$stages, xpd = TRUE, srt = 45, adj = 1)
    for(i in 1:d$nClusters){
        clusterProfilePlot$addLines(
            x = rep(d$stage_means[i], 2),
            y = ylim,
            col = d$colors[i],
            lwd = 1
        )
    }
    for(i in 1:d$nClusters){
        clusterProfilePlot$addLines(
            x = 1:d$nStages,
            y = d$clusters[i,],
            col = d$colors[i],
            lwd = 2
        )
    }
}
plotClusterHeatmap <- function(d) {
    par(mar = c(4, 4, 0, 0) + 0.1)
    clusterProfilePlot$initializeFrame(
        xlim = c(1, d$nStages),
        ylim = c(1, d$nClusters),
        xlab = "",
        ylab = "Cluster",
        xaxt = "n",
        yaxt = "n"
    )
    dt <- as.data.table(d$clusters)
    setnames(dt, d$stages)
    dt[, cluster := 1:nrow(dt)]
    dt <- melt(dt, id.vars = "cluster", measure.vars = d$stages, variable.name = "stage")
    dt$stage <- as.integer(dt$stage)
    dstr(dt)
    dt <- dt[, .(x = stage, y = cluster, z = value)]
    dstr(dt)
    mdiLevelPlot(
        dt,     # a data.table with at least columns x, y, and the column named by z.column
        xlim = c(0, d$nStages),   # the plot X-axis limits
        xinc = 1,   # the regular increment of the X-axis grid
        ylim = c(0, d$nClusters),   # the plot Y-axis limits
        yinc = 1,   # the regular increment of the Y-axis grid
        z.fn = function(z) z,   # function applied to z.column, per grid spot, to generate the output color
        z.column = "z", # the column in dt passed to z.fn, per grid spot
        settings = clusterProfilePlot$settings, 
        legendTitle = "RPKM"
    )
}
clusterProfilePlot <- staticPlotBoxServer(
    "clusterProfilePlot",
    margins = TRUE,
    settings = c(clusterSettings, mdiLevelPlotSettings),
    create = function() {
        d <- clusterAggregates()
        req(d)
        startSpinner(session, message = "rendering clusters")
        if(clusterProfilePlot$settings$get("Clustering", "Plots_Clusters_As") == "traces") {
            plotClusterTraces(d)
        } else {
            plotClusterHeatmap(d)
        }
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
