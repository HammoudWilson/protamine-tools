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
paRainbow <- function(n) rainbow(n / 0.9)[1:n] # eliminates magenta colors that are visually too close to red
titledMar <- c(4, 4, 1.7, 0.2) + 0.1
interactiveFrame <- list(
    Plot_Frame = list(
        Width_Inches = 4,
        Height_Inches = 3,
        Font_Size = 9
    )
)

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
stages <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTss_ab_initio(sourceId)$stages
})
observeEvent(sourceId(), {
    stages <- stages()
    req(stages)
    updateSelectInput(
        session, 
        "index_stage", 
        label = "Index Stage",
        choices = c("overlap_group", stages), 
        selected = "overlap_group"
    )
})
regions <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTss_dinuc_regions(sourceId, input$index_stage)
})

#----------------------------------------------------------------------
# interactive correlation plot
#----------------------------------------------------------------------
correlationPlotSettings <- list(
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
)
correlationPlotData <- reactive({
    stages <- stages()
    regions <- regions()
    req(stages, regions)
    startSpinner(session, message = "loading correlation")
    regions <- regions[order(stage_mean)]
    x_axis <- correlationPlotBox$settings$get("Correlation", "X_Axis")
    y_axis <- correlationPlotBox$settings$get("Correlation", "Y_Axis")
    list(
        x_axis = x_axis,
        y_axis = y_axis,
        xy = as.matrix(regions[, .SD, .SDcols = c(x_axis, y_axis)]),
        regions = regions,
        color = paRainbow(nrow(regions)), 
        stages = stages
    )
})
getCorrelationRange <- function(xy, col, stages) {
    if(col == "stage_mean") return(c(1, length(stages)))
    lim <- range(xy, na.rm = TRUE)
    lim[1] <- lim[1] - diff(lim) * 0.05 # extend range by 5% on each side
    lim[2] <- lim[2] + diff(lim) * 0.05
    lim
}
createCorrelationPlot <- function(settings, plot){
    d <- correlationPlotData()
    req(d)
    startSpinner(session, message = "rendering correlation")
    layout <- plot$initializePng(
        mar = titledMar
    ) %>% plot$initializeFrame(
        xlim = getCorrelationRange(d$xy[, 1], d$x_axis, d$stages),
        ylim = getCorrelationRange(d$xy[, 2], d$y_axis, d$stages),
        xlab = d$x_axis,
        ylab = d$y_axis,
        xaxs = "i",
        yaxs = "i",
        title = paste0(input$index_stage, " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    I <- sample(1:nrow(d$xy))
    plot$addPoints(
        x = d$xy[I, 1],
        y = d$xy[I, 2],
        col = d$color[I]
    )
    stopSpinner(session)
    plot$finishPng(layout)
}
correlationPlotBox <- mdiInteractivePlotBoxServer(
    "correlationPlotBox",
    click = TRUE,
    points = TRUE,
    settings = correlationPlotSettings,
    defaults = c(interactiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    )),
    create = function(...) createCorrelationPlot(..., correlationPlotBox)
)
observeEvent(correlationPlotBox$plot$click(), {
    setSelectedRegion(correlationPlotBox$plot$click()$coord, correlationPlotData())
})

#----------------------------------------------------------------------
# interactive umap plot
#----------------------------------------------------------------------
umapPlotData <- reactive({
    regions <- regions()
    req(regions)
    startSpinner(session, message = "loading umap data")
    umap1_col <- paste(input$rpkm_scaling, input$umap_metric, "umap1", sep = "_")
    umap2_col <- paste(input$rpkm_scaling, input$umap_metric, "umap2", sep = "_")
    regions <- regions[order(stage_mean)]
    list(
        xy      = regions[, .SD, .SDcols = c(umap1_col, umap2_col)],
        regions = regions,
        color   = paRainbow(nrow(regions))
    )
})
createUmapPlot <- function(settings, plot){
    d <- umapPlotData()
    req(d)
    startSpinner(session, message = "rendering umap plot")
    layout <- plot$initializePng(
        mar = titledMar
    ) %>% plot$initializeFrame(
        xlim = range(d$xy[, 1]),
        ylim = range(d$xy[, 2]),
        xlab = "UMAP 1",
        ylab = "UMAP 2",
        title = paste0(paste(input$index_stage, input$rpkm_scaling, input$umap_metric), " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    I <- sample(1:nrow(d$xy))
    plot$addPoints(
        x = d$xy[I][[1]],
        y = d$xy[I][[2]],
        col = d$color[I]
    )
    stopSpinner(session)
    plot$finishPng(layout)
}
umapPlotBox <- mdiInteractivePlotBoxServer(
    "umapPlotBox",
    click = TRUE,
    points = TRUE,
    defaults = c(interactiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.25
        )
    )),
    create = function(...) createUmapPlot(..., umapPlotBox) # a function or reactive that creates the plot as a png file using settings and helpers
)
observeEvent(umapPlotBox$plot$click(), {
    setSelectedRegion(umapPlotBox$plot$click()$coord, umapPlotData())
})

#----------------------------------------------------------------------
# profile plot by cluster/quantile
#----------------------------------------------------------------------
clusterProfilePlotSettings <- list(
    Clustering = list(
        Cluster_By = list(
            type = "selectInput",
            choices = c("quantile", "k_means"),
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
    regions <- regions()
    req(ai, regions)
    startSpinner(session, message = "aggregating clusters")
    stage_rpkm_cols   <- paste(ai$stages, "rpkm",   sep = '_')
    stage_scaled_cols <- paste(ai$stages, "scaled", sep = '_')
    stage_cols <- if(input$rpkm_scaling == "scaled") stage_scaled_cols else stage_rpkm_cols
    clusterCol <- if(clusterProfilePlot$settings$get("Clustering","Cluster_By") == "quantile") {
        "quantile"
    } else {
        paste(input$rpkm_scaling, "cluster", sep = "_")
    }
    clusters <- sapply(stage_cols,      function(col) regions[, mean(regions[[col]][.I]), keyby = clusterCol][[2]])
    weights  <- sapply(stage_rpkm_cols, function(col) regions[, mean(regions[[col]][.I]), keyby = clusterCol][[2]])
    clusterCounts <- regions[, .N, keyby = clusterCol]
    nClusters <- nrow(clusters)
    nStages   <- ncol(clusters)
    stage_means <- sapply(1:nClusters, function(i) matrixStats::weightedMedian(1:nStages, weights[i,]))
    if(clusterCol != "quantile") {
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
        colors = paRainbow(nClusters)
    )
})
addStageXAxis <- function(stages, ylim, side = 1) {
    axis(side = side, labels = FALSE)
    text(1:length(stages), ylim[1] - diff(ylim) / 8, stages, xpd = TRUE, srt = 45, adj = 1)
}
plotClusterTraces <- function(d) {
    ylim <- if(input$rpkm_scaling == "scaled") range(d$clusters) else c(0, max(d$clusters))
    ylim <- ylim * 1.05
    par(mar = titledMar)
    clusterProfilePlot$initializeFrame(
        xlim = c(1, d$nStages),
        ylim = ylim,
        xlab = "",
        ylab = if(input$rpkm_scaling == "scaled") "log2(stage RPKM / mean RPKM)" else "RPKM",
        xaxt = "n",
        title = paste0(input$index_stage, " " , input$rpkm_scaling, " (", format(nrow(d$clusters), big.mark = ","), ")")
    )
    addStageXAxis(d$ai$stages, ylim)
    if(input$rpkm_scaling == "scaled") abline(h = 0, col = CONSTANTS$plotlyColors$grey)
    for(i in 1:d$nClusters) abline(v = d$stage_means[i], col = d$colors[i])
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
    xlim <- c(0.5, d$nStages + 0.5)
    ylim <- c(0.5, d$nClusters + 0.5)
    par(mar = titledMar)
    clusterProfilePlot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = "Cluster",
        xaxt = "n",
        yaxt = "n",
        xaxs = "i",
        yaxs = "i",
        title = paste0(input$index_stage, " " , input$rpkm_scaling, " (", format(nrow(d$clusters), big.mark = ","), ")")
    )
    addStageXAxis(d$ai$stages, ylim)
    I <- d$nClusters:1 # invert so earliest is at the top
    dt <- as.data.table(d$clusters[I,])
    setnames(dt, d$ai$stages)
    dt[, cluster := 1:nrow(dt)]
    dt <- melt(
        dt, 
        id.vars = "cluster", 
        measure.vars = d$ai$stages, 
        variable.name = "stage"
    )
    dt$stage <- as.integer(dt$stage)
    dt <- dt[, .(x = stage, y = cluster, z = value)]
    mdiLevelPlot(
        dt,     # a data.table with at least columns x, y, and the column named by z.column
        xlim = xlim,   # the plot X-axis limits
        xinc = 1,   # the regular increment of the X-axis grid
        ylim = ylim,   # the plot Y-axis limits
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
    settings = c(clusterProfilePlotSettings, mdiLevelPlotSettings),
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
        dstr(d)
        startSpinner(session, message = "rendering profile")
        nStages <- length(ai$stages)
        stage_cols <- if(input$rpkm_scaling == "scaled") {
            paste(ai$stages, "scaled", sep = '_')
        } else {
            paste(ai$stages, "rpkm", sep = '_')
        }
        y <- unlist(d$region[, .SD, .SDcols = stage_cols])
        ylim <- if(input$rpkm_scaling == "scaled") range(y) else c(0, max(y))
        ylim <- ylim * 1.05
        par(mar = titledMar)
        regionProfilePlot$initializeFrame(
            xlim = c(1, nStages),
            ylim = ylim,
            xlab = "",
            ylab = "RPKM",
            xaxt = "n",
            title = paste0(d$region$chrom, ":", d$region$start0, "-", d$region$end1)
        )
        addStageXAxis(ai$stages, ylim)
        if(input$rpkm_scaling == "scaled") abline(h = 0, col = CONSTANTS$plotlyColors$grey)
        abline(v = d$region$stage_mean, col = d$color, lty = 1)
        regionProfilePlot$addLines(
            x = 1:nStages,
            y = y,
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
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        correlationPlotBox$settings$replace(bm$outcomes$correlationPlotBoxSettings)
        umapPlotBox$settings$replace(bm$outcomes$umapPlotBoxSettings)
        clusterProfilePlot$settings$replace(bm$outcomes$clusterProfilePlotSettings)
        regionProfilePlot$settings$replace(bm$outcomes$regionProfilePlotSettings)
        regionPlotBox$settings$replace(bm$outcomes$regionPlotBoxSettings)
    }
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
    outcomes = list(
        correlationPlotBoxSettings = correlationPlotBox$settings$all_,
        umapPlotBoxSettings = umapPlotBox$settings$all_,
        clusterProfilePlotSettings = clusterProfilePlot$settings$all_,
        regionProfilePlotSettings = regionProfilePlot$settings$all_,
        regionPlotBoxSettings = regionPlotBox$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
