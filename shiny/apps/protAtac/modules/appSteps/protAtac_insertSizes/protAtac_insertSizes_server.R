#----------------------------------------------------------------------
# server components for the protAtac_insertSizes appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
protAtac_insertSizesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'protAtac_insertSizes'
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
aggregateInsertSizes <- function(insertSizes, samples){
    stages <- unique(samples$stage)
    x <- do.call(cbind, lapply(stages, function(stage_){
        apply(insertSizes[, .SD, .SDcols = samples[stage == stage_, sample_name]], 1, sum, na.rm = TRUE)
    })) %>% as.data.table
    setnames(x, stages)
    x
}
normalizeInsertSizes <- function(this, ref){
    seriesNames <- colnames(this)
    x <- do.call(cbind, lapply(seriesNames, function(series){
        this[[series]] * (sum(this[[series]], na.rm = TRUE) / sum(ref[[series]], na.rm = TRUE))
    })) %>% as.data.table
    setnames(x, seriesNames)
    x
}
getInsertSizeData <- function(insertSizes, refType, samples, aggregate, normalize){
    this <- insertSizes[[refType]][, .SD, .SDcols = samples$sample_name]
    if(aggregate) this <- aggregateInsertSizes(this, samples)
    if(normalize){
        ref <-insertSizes$spike_in[, .SD, .SDcols = samples$sample_name]
        if(aggregate) ref <- aggregateInsertSizes(ref, samples)
        normalizeInsertSizes(this, ref)
    } else {
        this[, lapply(.SD, function(x) x / sum(x, na.rm = TRUE))]
    }
}
insertSizesPlot <- function(refType){
    plot <- staticPlotBoxServer(
        paste("insertSizesPlot", refType, sep = "_"),
        maxHeight = "400px",
        lines   = TRUE,
        legend  = TRUE,
        margins = TRUE,
        title   = TRUE,
        create = function() {
            sourceId <- sourceId()
            samples <- spermatidStages$selectedSamples()
            allSamples <- allSamples()
            req(sourceId, samples)
            isd <- paInsertSizes(sourceId)
            aggregate <- input$aggregateSamplesByStage
            normalize <- input$normalizeToSpikeIn
            if(refType == "spike_in" && normalize) req(FALSE)
            binSize <- isd$bin_size
            isd <- getInsertSizeData(isd$insertSizes, refType, samples, aggregate, normalize)
            seriesNames <- colnames(isd)
            colors <- if(aggregate) getStageColors(allSamples, samples) else getSampleColorsByStage(allSamples, samples)
            plot$initializeFrame(
                xlim = c(0, 700),
                ylim = c(0, max(isd)),
                xlab = "Insert Size (bp)",
                ylab = "Frequency"
            )
            for(series in seriesNames){
                plot$addLines( # addLines follows the same pattern, etc.
                    x = 1:binSize,
                    y = isd[[series]],
                    col = colors[series]
                )
            }
            abline(v = c(65, 125, 146), col = CONSTANTS$plotlyColors$grey) # intermediate and nucleosome size boundaries
            plot$addLegend(
                legend = seriesNames,
                col = colors,
                cex = 0.85
            )
            stopSpinner(session)
        }
    )
    plot
}
genomePlot  <- insertSizesPlot("genome")
spikeInPlot <- insertSizesPlot("spike_in")

# #----------------------------------------------------------------------
# # plot outputs
# #----------------------------------------------------------------------
# NRLLPlot <- staticPlotBoxServer(
#     "insertSizesPlot_NRLL",
#     maxHeight = "400px",
#     lines   = TRUE,
#     legend  = TRUE,
#     margins = TRUE,
#     title   = TRUE,
#     create = function() {
#         sourceId <- sourceId()
#         seg <- paSegmentation(sourceId)
#         samples <- spermatidStages$selectedSamples()
#         req(sourceId, seg, samples)
#         bd <- paBinData(sourceId)
#         d <- sapply(samples$sample_name, function(sample_name){
#             x <- seg$NRLL[[sample_name]][bd$bins$genome$excluded == 0]
#             x <- as.integer(round(x[x != 0.0] / 0.1)) * 0.1
#             data.table(x = x)[, .(y = .N), keyby = .(x)]
#         }, simplify = FALSE, USE.NAMES = TRUE)
#         colors <- sapply(samples$sample_name, function(sample_name_){
#             stageColors[samples[sample_name == sample_name_, stage]]
#         })
#         names(colors) <- samples$sample_name
#         NRLLPlot$initializeFrame(
#             xlim = c(-2, 1.5),
#             ylim = c(0, 0.25),
#             xlab = "Normalized Relative Log Likelihood (NRLL)",
#             ylab = "Frequency"
#         )
#         abline(v = 0, col = CONSTANTS$plotlyColors$grey)
#         for(sample_name in samples$sample_name){
#             NRLLPlot$addLines(
#                 x = d[[sample_name]]$x,
#                 y = d[[sample_name]]$y / sum(d[[sample_name]]$y, na.rm = TRUE),
#                 col = colors[sample_name]
#             )
#         }
#         NRLLPlot$addLegend(
#             legend = samples$sample_name,
#             col = colors,
#             cex = 0.8
#         )
#         stopSpinner(session)
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
        genomePlot$settings$replace(bm$outcomes$genomePlotSettings)
        spikeInPlot$settings$replace(bm$outcomes$spikeInPlotSettings)
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
        genomePlotSettings  = genomePlot$settings$all_,
        spikeInPlotSettings = spikeInPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
