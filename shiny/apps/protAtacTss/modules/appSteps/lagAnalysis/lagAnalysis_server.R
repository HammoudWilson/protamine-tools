#----------------------------------------------------------------------
# server components for the lagAnalysis appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
lagAnalysisServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'lagAnalysis'
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
lagData <- reactive({
    sourceId <- req(sourceId())
    file <- getSourceFilePath(sourceId, "lag_analysis")
    ld <- readRDS(file)
    list(
        xlim = range(ld$log10dist),
        xlab = "Inter-peak Distance (log10 bp)",
        ld = ld
    )
})

#----------------------------------------------------------------------
# variogram plot
#----------------------------------------------------------------------
createVariogramPlot <- function(settings, plotBox){
    ld <- lagData()
    req(ld)
    startSpinner(session, message = "rendering variogram")
    ldd <- ld$ld
    layout <- plotBox$initializePng(
        mar = titledMar
    ) %>% plotBox$initializeFrame(
        xlim = ld$xlim,
        ylim = c(0, max(ldd$var) * 1.05),
        xlab = ld$xlab,
        ylab = "Var(Stage Mean)",
        # xaxs = "i",
        # yaxs = "i",
        # title = paste0(chains$input$index_stage, " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    abline(h = 0, col = CONSTANTS$plotlyColors$black)
    variogramPlotBox$addPoints(
        x = ldd$log10dist, 
        y = ldd$var, 
        pch = 19, 
        col = CONSTANTS$plotlyColors$blue
    )
    stopSpinner(session)
    plotBox$finishPng(layout)
}
variogramPlotBox <- mdiInteractivePlotBoxServer(
    "variogramPlotBox",
    points = TRUE,
    defaults = c(paInteractiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    )),
    create = function(...) createVariogramPlot(..., variogramPlotBox)
)

#----------------------------------------------------------------------
# autocorrelation plot (sort of)
#----------------------------------------------------------------------
createAutocorrelationPlot <- function(settings, plotBox){
    ld <- lagData()
    req(ld)
    startSpinner(session, message = "rendering autocorrelation")
    ldd <- ld$ld
    layout <- plotBox$initializePng(
        mar = titledMar
    ) %>% plotBox$initializeFrame(
        xlim = ld$xlim,
        ylim = c(-1, 1),
        xlab = ld$xlab,
        ylab = "Mean Correlation",
        # xaxs = "i",
        # yaxs = "i",
        # title = paste0(chains$input$index_stage, " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    abline(h = 0, col = CONSTANTS$plotlyColors$black)
    segments(
        x0 = ldd$log10dist, 
        x1 = ldd$log10dist, 
        y0 = ldd$mean_corr - 1 * ldd$sd_corr, 
        y1 = ldd$mean_corr + 1 * ldd$sd_corr, 
        lwd = 0.5, 
        col = CONSTANTS$plotlyColors$grey
    )
    autocorrelationPlotBox$addPoints(
        x = ldd$log10dist, 
        y = ldd$mean_corr, 
        pch = 19, 
        col = CONSTANTS$plotlyColors$blue
    )
    stopSpinner(session)
    plotBox$finishPng(layout)
}
autocorrelationPlotBox <- mdiInteractivePlotBoxServer(
    "autocorrelationPlotBox",
    points = TRUE,
    defaults = c(paInteractiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    )),
    create = function(...) createAutocorrelationPlot(..., autocorrelationPlotBox)
)

#----------------------------------------------------------------------
# autocorrelation plot (sort of)
#----------------------------------------------------------------------
createDiffPlot <- function(settings, plotBox){
    ld <- lagData()
    req(ld)
    startSpinner(session, message = "rendering stage mean diff")
    ldd <- ld$ld
    layout <- plotBox$initializePng(
        mar = titledMar
    ) %>% plotBox$initializeFrame(
        xlim = ld$xlim,
        ylim = c(0, 1.5),
        xlab = ld$xlab,
        ylab = "Mean Stage Mean Delta",
        # xaxs = "i",
        # yaxs = "i",
        # title = paste0(chains$input$index_stage, " (", format(nrow(d$xy), big.mark = ","), ")"),
        cex.main = 1.05
    )
    abline(h = 0, col = CONSTANTS$plotlyColors$black)
    segments(
        x0 = ldd$log10dist, 
        x1 = ldd$log10dist, 
        y0 = ldd$mean_diff - 1 * ldd$sd_diff, 
        y1 = ldd$mean_diff + 1 * ldd$sd_diff, 
        lwd = 0.5, 
        col = CONSTANTS$plotlyColors$grey
    )
    diffPlotBox$addPoints(
        x = ldd$log10dist, 
        y = ldd$mean_diff, 
        pch = 19, 
        col = CONSTANTS$plotlyColors$blue
    )
    stopSpinner(session)
    plotBox$finishPng(layout)
}
diffPlotBox <- mdiInteractivePlotBoxServer(
    "diffPlotBox",
    points = TRUE,
    defaults = c(paInteractiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    )),
    create = function(...) createDiffPlot(..., diffPlotBox)
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        variogramPlotBox$settings$replace(bm$outcomes$variogramPlotBoxSettings)
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
        variogramPlotBoxSettings = variogramPlotBox$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
