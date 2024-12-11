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

# #----------------------------------------------------------------------
# # plot outputs
# #----------------------------------------------------------------------
# insertSizesPlot <- function(refType){
#     plot <- staticPlotBoxServer(
#         paste("insertSizesPlot", refType, sep = "_"),
#         maxHeight = "400px",
#         lines   = TRUE,
#         legend  = TRUE,
#         margins = TRUE,
#         title   = TRUE,
#         create = function() {
#             sourceId <- sourceId()
#             samples <- spermatidStages$selectedSamples()
#             req(sourceId, samples)
#             isd <- paInsertSizes(sourceId)
#             aggregate <- settings$get("Insert_Sizes","Aggregate_Samples_By_Stage")
#             normalize <- settings$get("Insert_Sizes","Normalize_To_Spike_In")
#             if(refType == "spike_in" && normalize) req(FALSE)
#             binSize <- isd$bin_size
#             isd <- getInsertSizeData(isd$insertSizes, refType, samples, aggregate, normalize)
#             seriesNames <- colnames(isd)
#             colors <- stageColors[if(aggregate) seriesNames else samples$stage]
#             # colors <- 1:length(seriesNames)
#             names(colors) <- seriesNames
#             plot$initializeFrame(
#                 xlim = c(0, 700),
#                 ylim = c(0, max(isd)),
#                 xlab = "Insert Size (bp)",
#                 ylab = "Frequency"
#             )
#             for(series in seriesNames){
#                 plot$addLines( # addLines follows the same pattern, etc.
#                     x = 1:binSize,
#                     y = isd[[series]],
#                     col = colors[series]
#                 )
#             }
#             abline(v = 146, col = CONSTANTS$plotlyColors$grey) # nucleosome size
#             plot$addLegend(
#                 legend = seriesNames,
#                 col = colors,
#                 cex = 0.85
#             )
#             stopSpinner(session)
#         }
#     )
#     plot
# }
# genomePlot  <- insertSizesPlot("genome")
# spikeInPlot <- insertSizesPlot("spike_in")

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
plotDistributions <- function(plot, scoreType, scores, colors, message){
    startSpinner(session, message = message)
    maxY <- max(unlist(sapply(names(scores), function(seriesName) scores[[seriesName]]$dist$y)), na.rm = TRUE)
    plot$initializeFrame(
        xlim = scoreType$lim,
        ylim = c(0, maxY * 1.05),
        xlab = scoreType$name,
        ylab = "Frequency"
    )
    abline(v = 0, col = CONSTANTS$plotlyColors$grey)
    for(seriesName in names(scores)){
        plot$addLines(
            x = scores[[seriesName]]$dist,
            col = colors[seriesName]
        )
        abline(v = scores[[seriesName]]$mean, col = colors[seriesName])
    }
    plot$addLegend(
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

        dstr(getStageTypeScores(sourceId, input$scoreType, samples))
        dstr(getStageTypeColors(sourceId, allSamples, samples))
        # req(FALSE)

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
        # scoreLevel <- getScoreLevel(input$scoreType)
        sourceId <- sourceId()
        req(sourceId)

        dstr(getStageTypeDeltaScores(sourceId, input$scoreType))
        dstr(list(stageType_delta = CONSTANTS$plotlyColors$black))
        # req(FALSE)

        plotDistributions(
            plot        = deltaDistributionPlot, 
            scoreType   = getScoreType(sourceId, input$scoreType), 
            scores      = getStageTypeDeltaScores(sourceId, input$scoreType),
            colors      = c(stageType_delta = CONSTANTS$plotlyColors$black),
            message     = "plotting delta"
        )
    }
)
# deltaZDistributionPlot <- staticPlotBoxServer(
#     "deltaZDistributionPlot",
#     maxHeight = "400px",
#     lines   = TRUE,
#     legend  = TRUE,
#     margins = TRUE,
#     title   = TRUE,
#     create = function() {
#         # scoreLevel <- getScoreLevel(input$scoreType)
#         sourceId <- sourceId()
#         req(sourceId)

#         dstr(getStageTypeDeltaScores(sourceId, input$scoreType))
#         dstr(list(stageType_delta = CONSTANTS$plotlyColors$black))
#         # req(FALSE)

#         plotDistributions(
#             plot        = deltaZDistributionPlot, 
#             scoreType   = getScoreType(sourceId, input$scoreType), 
#             scores      = getStageTypeDeltaScores(sourceId, input$scoreType),
#             colors      = c(stageType_delta = CONSTANTS$plotlyColors$black),
#             message     = "plotting deltaZ"
#         )
#     }
# )

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
