#----------------------------------------------------------------------
# server components for the dinucChains appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
dinucChainsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'dinucChains'
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
chains   <- dinucChainsSelectorBoxServer("dinucChains", sourceId)
clusterProfilePlotBox <- clusterProfilePlotBoxServer("clusterProfilePlotBox", chains)
intervalExpansion     <- intervalExpansionServer("intervalExpansion", chains)

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
    intervals <- chains$intervals()
    stages    <- chains$stages()
    req(intervals, stages)
    startSpinner(session, message = "loading correlation")
    intervals <- intervals[order(stage_mean)]
    x_axis <- correlationPlotBox$settings$get("Correlation", "X_Axis")
    y_axis <- correlationPlotBox$settings$get("Correlation", "Y_Axis")
    list(
        x_axis = x_axis,
        y_axis = y_axis,
        xy = as.matrix(intervals[, .SD, .SDcols = c(x_axis, y_axis)]),
        intervals = intervals,
        color = paRainbow(nrow(intervals)),
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
createCorrelationPlot <- function(settings, plotBox){
    d <- correlationPlotData()
    req(d)
    startSpinner(session, message = "rendering correlation")
    layout <- plotBox$initializePng(
        mar = titledMar
    ) %>% plotBox$initializeFrame(
        xlim = getCorrelationRange(d$xy[, 1], d$x_axis, d$stages),
        ylim = getCorrelationRange(d$xy[, 2], d$y_axis, d$stages),
        xlab = d$x_axis,
        ylab = d$y_axis,
        xaxs = "i",
        yaxs = "i",
        title = paste0(chains$input$index_stage, " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    randomOrder <- sample(1:nrow(d$xy))
    plotBox$addPoints(
        x = d$xy[randomOrder, 1],
        y = d$xy[randomOrder, 2],
        col = d$color[randomOrder]
    )
    stopSpinner(session)
    plotBox$finishPng(layout)
}
correlationPlotBox <- mdiInteractivePlotBoxServer(
    "correlationPlotBox",
    click = TRUE,
    points = TRUE,
    settings = correlationPlotSettings,
    defaults = c(paInteractiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    )),
    create = function(...) createCorrelationPlot(..., correlationPlotBox)
)
observeEvent(correlationPlotBox$plot$click(), {
    intervalExpansion$setSelectedInterval(correlationPlotBox$plot$click()$coord, correlationPlotData())
})

#----------------------------------------------------------------------
# interactive umap plot
#----------------------------------------------------------------------
umapPlotData <- reactive({
    intervals <- chains$passedIntervals()
    req(intervals)
    startSpinner(session, message = "loading umap data")
    umap1_col <- paste(chains$input$rpkm_scaling, chains$input$umap_metric, "umap1", sep = "_")
    umap2_col <- paste(chains$input$rpkm_scaling, chains$input$umap_metric, "umap2", sep = "_")
    intervals <- intervals[order(stage_mean)]
    list(
        xy        = as.matrix(intervals[, .SD, .SDcols = c(umap1_col, umap2_col)]),
        intervals = intervals,
        color     = paRainbow(nrow(intervals))
    )
})
createUmapPlot <- function(settings, plotBox){
    d <- umapPlotData()
    req(d)
    startSpinner(session, message = "rendering umap plot")
    layout <- plotBox$initializePng(
        mar = titledMar
    ) %>% plotBox$initializeFrame(
        xlim = range(d$xy[, 1]),
        ylim = range(d$xy[, 2]),
        xlab = "UMAP 1",
        ylab = "UMAP 2",
        title = paste0(paste(
            chains$input$index_stage, 
            chains$input$rpkm_scaling, 
            chains$input$umap_metric
        ), " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    randomOrder <- sample(1:nrow(d$xy))
    plotBox$addPoints(
        x = d$xy[randomOrder, 1],
        y = d$xy[randomOrder, 2],
        col = d$color[randomOrder]
    )
    stopSpinner(session)
    plotBox$finishPng(layout)
}
umapPlotBox <- mdiInteractivePlotBoxServer(
    "umapPlotBox",
    click = TRUE,
    points = TRUE,
    defaults = c(paInteractiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.25
        )
    )),
    create = function(...) createUmapPlot(..., umapPlotBox) # a function or reactive that creates the plot as a png file using settings and helpers
)
observeEvent(umapPlotBox$plot$click(), {
    intervalExpansion$setSelectedInterval(umapPlotBox$plot$click()$coord, umapPlotData())
})

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
        clusterProfilePlotBox$settings$replace(bm$outcomes$clusterProfilePlotBoxSettings)
        intervalExpansion$profilePlotBox$settings$replace(bm$outcomes$profilePlotBoxSettings)
        intervalExpansion$browserPlotBox$settings$replace(bm$outcomes$browserPlotBoxSettings)
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
        clusterProfilePlotBoxSettings = clusterProfilePlotBox$settings$all_,
        profilePlotBoxSettings = intervalExpansion$profilePlotBox$settings$all_,
        browserPlotBoxSettings = intervalExpansion$browserPlotBox$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
