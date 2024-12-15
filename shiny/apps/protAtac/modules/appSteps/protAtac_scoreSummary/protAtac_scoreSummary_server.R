#----------------------------------------------------------------------
# server components for the protAtac_scoreSummary appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
protAtac_scoreSummaryServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'protAtac_scoreSummary'
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
# data sources
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "dataSourceTable", 
    selection = "single"
) 
spermatidStages <- spermatidStageTableServer(
    "spermatidStageTable", 
    sourceId, # a reactive that returns the id of one selected source
    selection = "multiple"
)
allSamples <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    paBinData(sourceId)$samples
})

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
plotDistributions <- function(plot, scoreType, scores, colors, message, 
                              legend = TRUE, isDelta = FALSE){
    startSpinner(session, message = message)
    maxY <- max(unlist(sapply(names(scores), function(seriesName) scores[[seriesName]]$dist$y)), na.rm = TRUE)
    label <- paste(scoreType$label, scoreType$unit)
    plot$initializeFrame(
        xlim = scoreType[[if(isDelta) "deltaLim" else "valueLim"]],
        ylim = c(0, maxY * 1.05),
        xlab = if(isDelta) paste(label, "Delta") else label,
        ylab = "Frequency"
    )
    abline(v = 0, col = CONSTANTS$plotlyColors$grey)
    for(seriesName in names(scores)){
        plot$addLines(
            x = scores[[seriesName]]$dist,
            col = colors[seriesName]
        )
        abline(v = scores[[seriesName]]$median, col = colors[seriesName])
    }
    if(legend) plot$addLegend(
        legend = names(scores),
        col = colors,
        cex = 0.8
    )
    stopSpinner(session)
}
sampleDistributionPlot <- staticPlotBoxServer(
    "sampleDistributionPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        scoreLevel <- getScoreLevel(input$scoreType)
        if(scoreLevel != 'sample') req(FALSE)
        sourceId <- sourceId()
        req(sourceId)
        samples <- spermatidStages$selectedSamples()
        allSamples <- allSamples()
        plotDistributions(
            plot        = sampleDistributionPlot, 
            scoreType   = getScoreType(sourceId, input$scoreType), 
            scores      = getSampleScores(sourceId, input$scoreType, samples),
            colors      = getSampleColorsByStage(allSamples, samples),
            message     = "plotting samples"
        )
    }
)
stageDistributionPlot <- staticPlotBoxServer(
    "stageDistributionPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        scoreLevel <- getScoreLevel(input$scoreType)
        if(scoreLevel != 'sample') req(FALSE)
        sourceId <- sourceId()
        req(sourceId)
        samples <- spermatidStages$selectedSamples()
        allSamples <- allSamples()
        plotDistributions(
            plot        = stageDistributionPlot, 
            scoreType   = getScoreType(sourceId, input$scoreType), 
            scores      = getStageScores(sourceId, input$scoreType, samples),
            colors      = getStageColors(allSamples, samples),
            message     = "plotting stages"
        )
    }
)
stageTypeDistributionPlot <- staticPlotBoxServer(
    "stageTypeDistributionPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        scoreLevel <- getScoreLevel(input$scoreType)
        if(scoreLevel != 'sample') req(FALSE)
        sourceId <- sourceId()
        req(sourceId)
        samples <- spermatidStages$selectedSamples()
        allSamples <- allSamples()
        plotDistributions(
            plot        = stageTypeDistributionPlot, 
            scoreType   = getScoreType(sourceId, input$scoreType), 
            scores      = getStageTypeScores(sourceId, input$scoreType, samples),
            colors      = getStageTypeColors(sourceId, allSamples, samples),
            message     = "plotting stageTypes"
        )
    }
)
deltaDistributionPlot <- staticPlotBoxServer(
    "deltaDistributionPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        scoreLevel <- getScoreLevel(input$scoreType)
        sourceId <- sourceId()
        req(sourceId)
        plotDistributions(
            plot        = deltaDistributionPlot, 
            scoreType   = getScoreType(sourceId, input$scoreType), 
            scores      = getStageTypeDeltaScores(sourceId, input$scoreType, clean = TRUE),
            colors      = c(stageType_delta = CONSTANTS$plotlyColors$black),
            message     = "plotting delta",
            legend      = FALSE,
            isDelta     = scoreLevel == "sample"
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
        # spikeInPlot$settings$replace(bm$outcomes$spikeInPlotSettings)
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
        # spikeInPlotSettings = spikeInPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
