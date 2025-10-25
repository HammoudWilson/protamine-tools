#----------------------------------------------------------------------
# server components for the scoreCorrelation appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
scoreCorrelationServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'scoreCorrelation'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data sources
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "dataSourceTable", 
    selection = "single"
)

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
randomBinI <- reactive({
    sourceId <- req(sourceId())
    req(sourceId)
    bins <- paScores_bins(sourceId)
    includedAutosomeBins <- getIncludedAutosomeBins_scores(bins)
    binI <- paScores_bins(sourceId)[, which(includedAutosomeBins & !is.na(stgm))]
    sample(binI, settings$get("Score_Correlation", "Random_Bin_Count"))
})
getAxisScoreType <- function(sourceId, scoreTypeName, isZAxis = FALSE){
    scoreLevel <- getScoreLevel(scoreTypeName)
    randomBinI <- req(randomBinI())
    # metadata <- getStageTypeDeltaMetadata(sourceId, scoreTypeName, cleanDist = TRUE)
    scores <- if(scoreLevel == "sample") {
        getStageTypeDelta_allBins(sourceId, scoreTypeName)[randomBinI]
    } else {
        getGenomeScores(sourceId, scoreTypeName, randomBinI, stgmQuantile = isZAxis, hicQuantile = isZAxis) # do not apply quantile transform to genome scores here
    }
    list(
        scoreTypeName = scoreTypeName, 
        scoreType     = getScoreType(scoreTypeName),
        scoreLevel    = scoreLevel,
        isDelta       = scoreLevel == "sample",
        score         = scores
    )
}
xAxis <- reactive({
    sourceId <- req(sourceId())
    getAxisScoreType(sourceId, input$xScoreType)
})
yAxis <- reactive({
    sourceId <- req(sourceId())
    req(sourceId)
    getAxisScoreType(sourceId, input$yScoreType)
})
zAxis <- reactive({
    sourceId <- req(sourceId())
    req(sourceId)
    getAxisScoreType(sourceId, input$zScoreType, TRUE)
})
plotCorrelation <- function(plot, xAxis, yAxis, zAxis){
    startSpinner(session, message = "plotting correlation")
    xLabel <- paste(xAxis$scoreType$label, xAxis$scoreType$unit)
    yLabel <- paste(yAxis$scoreType$label, yAxis$scoreType$unit)
    plot$initializeFrame(
        xlim = xAxis$scoreType[[if(xAxis$isDelta) "deltaLim" else "valueLim"]],
        ylim = yAxis$scoreType[[if(yAxis$isDelta) "deltaLim" else "valueLim"]],
        xlab = if(xAxis$isDelta) paste(xLabel, "(Round - Elong)") else xLabel,
        ylab = if(yAxis$isDelta) paste(yLabel, "(Round - Elong)") else yLabel
    )
    plot$addPoints(
        x = xAxis$score,
        y = yAxis$score,
        # col = CONSTANTS$plotlyColors$blue %>% addAlphaToColor(0.1),
        col = getSeriesSummaryColors(zAxis$scoreTypeName, zAxis$score, 
                                     list(Max_Z_Score = 2, Min_Txn_Log10_CPM = -3, Max_Txn_Log10_CPM = 3)),
        pch = 16,
        cex = 0.5
    )
    # overplot a linear model fit
    fit <- lm(y ~ x, data.frame(x = xAxis$score, y = yAxis$score))
    xFit <- seq(    
        from = min(xAxis$score, na.rm = TRUE), 
        to   = max(xAxis$score, na.rm = TRUE), 
        length.out = 100
    )
    yFit <- predict(fit, newdata = data.frame(x = xFit))
    plot$addLines(
        x = xFit,
        y = yFit,
        col = CONSTANTS$plotlyColors$red,
        lwd = 1
    )
    stopSpinner(session)
}
correlationPlot <- staticPlotBoxServer(
    "correlationPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        plotCorrelation(
            plot     = correlationPlot, 
            xAxis    = xAxis(),
            yAxis    = yAxis(),
            zAxis    = zAxis()
        )
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        # genomePlot$settings$replace(bm$outcomes$genomePlotSettings)
        correlationPlot$settings$replace(bm$outcomes$correlationPlotSettings)
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
        # genomePlotSettings  = genomePlot$settings$all_,
        correlationPlotSettings = correlationPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
