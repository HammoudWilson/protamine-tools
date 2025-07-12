#----------------------------------------------------------------------
# server components for the clusterProfilePlotBox widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
clusterProfilePlotBoxServer <- function(id, sourceId, dinuc) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'clusterProfilePlotBox'
# settings <- activateMdiHeaderLinks( # uncomment as needed
#     session,
#     url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
#     dir = getAppStepDir(module), # for terminal emulator
#     envir = environment(), # for R console
#     baseDirs = getAppStepDir(module), # for code viewer/editor
#     settings = id, # for step-level settings
#     immediate = TRUE # plus any other arguments passed to settingsServer()
# )

#----------------------------------------------------------------------
# cluster profile plot
#----------------------------------------------------------------------
clusterProfilePlotSettings <- list(
    Clustering = list(
        # Cluster_By = list(
        #     type = "selectInput",
        #     choices = c("quantile", "k_means"),
        #     value = "quantile"
        # ),
        Plots_Clusters_As = list(
            type = "selectInput",
            choices = c("traces", "heatmap"),
            value = "traces"
        )
    )
)
stageColMean <- function(d){ # drop region+stage where RPKM==0 led to scaled log2 fold-change of -Inf, ~always late_ES
    if(dinuc$input$rpkm_scaling == "scaled") d <- d[is.finite(d)]
    if(length(d) == 0) NA_real_ else mean(d, na.rm = TRUE)
}
clusterAggregates <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    ai <- paTss_ab_initio(sourceId)
    req(ai)
    Cluster_By <- "quantile" # plot$settings$get("Clustering", "Cluster_By")
    regions <- if(Cluster_By == "quantile") {
        if(dinuc$input$include_unpassed_regions){
            clusterCol <- "quantile_unfiltered"
            dinuc$regions()
        } else {
            clusterCol <- "quantile_filtered"
            dinuc$passedRegions()
        }
    } else {
        clusterCol <- paste(dinuc$input$rpkm_scaling, "cluster", sep = "_")
        dinuc$passedRegions() # kmeans clustering was not performed with unpassed regions
    }
    req(regions)
    startSpinner(session, message = "aggregating clusters")
    stage_rpkm_cols   <- paste(ai$stages, "rpkm",   sep = '_')
    stage_scaled_cols <- paste(ai$stages, "scaled", sep = '_')
    stage_cols <- if(dinuc$input$rpkm_scaling == "scaled") stage_scaled_cols else stage_rpkm_cols
    clusters <- sapply(stage_cols,      function(col) regions[, stageColMean(regions[[col]][.I]), keyby = clusterCol][[2]])
    rpkms    <- sapply(stage_rpkm_cols, function(col) regions[, stageColMean(regions[[col]][.I]), keyby = clusterCol][[2]])
    clusterCounts <- regions[, .N, keyby = clusterCol]
    nClusters <- nrow(clusters)
    nStages   <- ncol(clusters)
    stage_means <- sapply(1:nClusters, function(i) matrixStats::weightedMedian(1:nStages, rpkms[i,], interpolate = TRUE, na.rm = TRUE))
    I <- order(stage_means)
    clusters <- clusters[I, ]
    stage_means <- stage_means[I]
    clusterCounts <- clusterCounts[I]
    list(
        ai = ai,
        clusters = clusters,
        clusterCounts = clusterCounts,
        stage_means = stage_means,
        nClusters = nClusters,
        nStages = nStages,
        colors = paRainbow(nClusters),
        regions = regions,
        clusterCol = clusterCol,
        clusterType = if(Cluster_By == "quantile") "Quantile (5%)" else "K-means Cluster"
    )
})
plotClusterTraces <- function(d) {
    ylim <- if(dinuc$input$rpkm_scaling == "scaled") range(d$clusters, na.rm = TRUE) 
            else c(0, max(d$clusters, na.rm = TRUE))
    ylim <- ylim * 1.05
    par(mar = titledMar)
    plot$initializeFrame(
        xlim = c(1, d$nStages),
        ylim = ylim,
        xlab = "",
        ylab = if(dinuc$input$rpkm_scaling == "scaled") "log2(stage RPKM / mean RPKM)" else "RPKM",
        xaxt = "n",
        title = paste0(dinuc$input$index_stage, " " , dinuc$input$rpkm_scaling, " (", format(sum(d$clusterCounts), big.mark = ","), ")")
    )
    addStageXAxis(d$ai$stages, ylim)
    if(dinuc$input$rpkm_scaling == "scaled") abline(h = 0, col = CONSTANTS$plotlyColors$grey)
    for(i in 1:d$nClusters) abline(v = d$stage_means[i], col = d$colors[i])
    for(i in 1:d$nClusters){
        plot$addLines(
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
    plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = "Cluster",
        xaxt = "n",
        yaxt = "n",
        xaxs = "i",
        yaxs = "i",
        title = paste0(dinuc$input$index_stage, " " , dinuc$input$rpkm_scaling, " (", format(nrow(d$clusters), big.mark = ","), ")")
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
        settings = plot$settings, 
        legendTitle = "RPKM"
    )
}
plot <- staticPlotBoxServer(
    "plot",
    margins = TRUE,
    settings = c(clusterProfilePlotSettings, mdiLevelPlotSettings),
    create = function() {
        d <- clusterAggregates()
        req(d)
        startSpinner(session, message = "rendering clusters")
        if(plot$settings$get("Clustering", "Plots_Clusters_As") == "traces") {
            plotClusterTraces(d)
        } else {
            plotClusterHeatmap(d)
        }
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
plot$data <- clusterAggregates
plot

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
