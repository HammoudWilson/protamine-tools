#----------------------------------------------------------------------
# server components for the cuttagSummary appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
cuttagSummaryServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'cuttagSummary'
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
cuttagSamplesTable <- cuttagSamplesTableServer(
    "cuttagSamplesTable", 
    samples = paCutTag_samples_reactive(sourceId),
    selection = "multiple"
)
selectedSamples <- reactive({
    cuttagSamplesTable$selectedSamples()[antibody_target == input$scoreType]
})
selectedSampleNames <- reactive({
    selectedSamples()$sample_name
})
allSamples <- reactive({
    cuttagSamplesTable$allSamples()
})
env <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paCollate_env(sourceId)
})

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
plotDistributions <- function(plot, scoreType, metadata, colors, message, 
                              samples, legend = TRUE){
    startSpinner(session, message = message)
    maxY <- max(unlist(sapply(names(metadata), function(seriesName) metadata[[seriesName]]$dist$y)), na.rm = TRUE)
    label <- paste(scoreType$label, scoreType$unit)
    plot$initializeFrame(
        xlim = scoreType$valueLim,
        # ylim = log10(c(1e-4, maxY * 1.05)), # c(0, maxY * 1.05),
        ylim = c(0, maxY * 1.05),
        xlab = label,
        ylab = "Frequency"
    )
    abline(v = 0, col = CONSTANTS$plotlyColors$grey)
    for(seriesName in names(metadata)){
        if(!is.null(samples)){
            sample <- samples[sample_name == seriesName]
            x <- 0:scoreType$valueLim[2]
            plot$addLines(
                x = x,
                # y = log10(dpois(x, lambda = sample$reads_per_bin)),
                y = dpois(x, lambda = sample$reads_per_bin),
                col = colors[seriesName]
            )
            abline(v = sample$reads_per_bin, col = colors[seriesName], lty = 2)
        }
        plot$addPoints(
            x = metadata[[seriesName]]$dist$x,
            # y = log10(metadata[[seriesName]]$dist$y),
            y = metadata[[seriesName]]$dist$y,
            col = colors[seriesName]
        )
        abline(v = metadata[[seriesName]]$median, col = colors[seriesName])
    }
    if(legend) plot$addLegend(
        legend = names(metadata),
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
        sourceId <- sourceId()
        req(sourceId)
        samples <- selectedSamples()
        allSamples <- allSamples()
        plotDistributions(
            plot        = sampleDistributionPlot, 
            scoreType   = getScoreType(input$scoreType), 
            metadata    = getCutTagSampleMetadata(sourceId, input$scoreType, samples),
            colors      = getSampleColorsByStage(allSamples, samples),
            message     = "plotting samples",
            samples     = samples
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
        sourceId <- sourceId()
        req(sourceId)
        samples <- selectedSamples()
        allSamples <- allSamples()
        plotDistributions(
            plot        = stageDistributionPlot, 
            scoreType   = getScoreType(input$scoreType), 
            metadata    = getCutTagStageMetadata(sourceId, input$scoreType, samples),
            colors      = getStageColors(allSamples, samples),
            message     = "plotting stages",
            samples     = NULL
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
        sampleDistributionPlot$settings$replace(bm$outcomes$sampleDistributionPlotSettings)
        stageDistributionPlot$settings$replace(bm$outcomes$stageDistributionPlotSettings)
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
        sampleDistributionPlotSettings  = sampleDistributionPlot$settings$all_,
        stageDistributionPlotSettings = stageDistributionPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
