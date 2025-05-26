#----------------------------------------------------------------------
# server components for the insertSizes appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
insertSizesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'insertSizes'
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
spermatidSamplesTable <- spermatidSamplesTableServer(
    "spermatidSamplesTable", 
    sourceId,
    selection = "multiple"
)
selectedSamples <- reactive({
    spermatidSamplesTable$selectedSamples()
})
allSamples <- reactive({
    spermatidSamplesTable$allSamples()
})

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
aggregateInsertSizes <- function(isdt, samples){
    # sum insert size counts for all samples in the same stage
    stages <- unique(samples$stage)
    x <- do.call(cbind, lapply(stages, function(stage_){
        apply(isdt[, .SD, .SDcols = samples[stage == stage_, sample_name]], 1, sum, na.rm = TRUE)
    })) %>% as.data.table
    setnames(x, stages)
    x
}
normalizeInsertSizes <- function(isdt_this, isdt_ref){
    seriesNames <- colnames(isdt_this)
    x <- do.call(cbind, lapply(seriesNames, function(series){
        nInserts_this <- sum(isdt_this[[series]], na.rm = TRUE)
        nInserts_ref  <- sum(isdt_ref[[series]], na.rm = TRUE)
        dmsg(nInserts_this)
        dmsg(nInserts_ref)
        isdt_this[[series]] * (nInserts_this / nInserts_ref)
    })) %>% as.data.table
    setnames(x, seriesNames)
    x
}
getInsertSize_refType <- function(insertSizes, refType, samples){
    # return a data.table of all insert sizes from 1 to MAX_INSERT_SIZE, one column per requested sample
    sapply(samples$sample_name, function(sample_name){
        c(
            rep(0, MIN_INSERT_SIZE - 1),
            rowSums(insertSizes[[sample_name]][[refType]])
        )
    }, simplify = FALSE, USE.NAMES = TRUE) %>% as.data.table
}
getInsertSizeData <- function(insertSizes, refType, samples, aggregate, normalize){
    isdt_this <- getInsertSize_refType(insertSizes, refType, samples)
    if(aggregate) isdt_this <- aggregateInsertSizes(isdt_this, samples)
    if(normalize){
        isdt_ref <-getInsertSize_refType(insertSizes, "spike_in", samples)
        if(aggregate) isdt_ref <- aggregateInsertSizes(isdt_ref, samples)
        normalizeInsertSizes(isdt_this, isdt_ref)
    } else {
        isdt_this[, lapply(.SD, function(x) x / sum(x, na.rm = TRUE))]
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
            selectedSamples <- selectedSamples()
            allSamples <- allSamples()
            req(sourceId, selectedSamples)
            isd <- paInsertSizes(sourceId)
            aggregate <- input$aggregateByStage
            normalize <- input$normalizeToSpikeIn
            if(refType == "spike_in" && normalize) req(FALSE)
            if(refType == "combined" && nrow(selectedSamples) > 1) req(FALSE)
            isd <- getInsertSizeData(isd$insertSizes, refType, selectedSamples, aggregate, normalize)
            seriesNames <- colnames(isd)
            colors <- if(aggregate) getStageColors(allSamples, selectedSamples) 
                      else if(length(unique(selectedSamples$stage)) == 1) getSampleColorsBySample(selectedSamples)
                      else getSampleColorsByStage(allSamples, selectedSamples)
            plot$initializeFrame(
                xlim = c(MIN_INSERT_SIZE, MAX_INSERT_SIZE),
                ylim = c(0, max(isd, na.rm = TRUE)),
                xlab = "Insert Size (bp)",
                ylab = "Frequency",
                xaxs = "i"
            )
            # abline(v = seq(100, 600, 100), col = "grey80") # nucleosome size boundaries
            abline(v = MAPPABILITY_KMER_LENGTHS - 0.5, col = "grey80")
            for(series in seriesNames){
                plot$addLines( # addLines follows the same pattern, etc.
                    x = MIN_INSERT_SIZE:MAX_INSERT_SIZE,
                    y = isd[[series]][MIN_INSERT_SIZE:MAX_INSERT_SIZE],
                    col = colors[series]
                )
            }
            plot$addLegend(
                legend = names(colors),
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
# spikeInPlot <- insertSizesPlot("combined")

#----------------------------------------------------------------------
# insert size vs. GC matrix heatmap plots
#----------------------------------------------------------------------
getInsertSizeMatrix <- function(insertSizes, refType, sampleNames){
    ism <- insertSizes[[sampleNames[1]]][[refType]]
    if(length(sampleNames) > 1){
        for (i in 2:length(sampleNames)){
            ism <- ism + insertSizes[[sampleNames[i]]][[refType]]
        }
    }
    ism
}
insertSizeMatrixPlot <- function(refType){
    plot <- staticPlotBoxServer(
        paste("insertSizeMatrix", refType, sep = "_"),
        maxHeight = "400px",
        create = function() {
            sourceId <- sourceId()
            selectedSamples <- selectedSamples()
            allSamples <- allSamples()
            req(sourceId, selectedSamples)
            isd <- paInsertSizes(sourceId)
            ism <- getInsertSizeMatrix(isd$insertSizes, refType, selectedSamples$sample_name)
            ism <- if(input$scaleHeatmapByInsertSize){
                apply(ism, 1, function(x) x / sum(x, na.rm = TRUE))
            } else {
                t(ism)
            }
            gcPercents <- 25:75
            minInsertSize <- if(input$nucleosomesOnly) 120 else MIN_INSERT_SIZE
            insertSizes <- minInsertSize:MAX_INSERT_SIZE
            nDummySizes <- MIN_INSERT_SIZE - 1
            ism <- cbind(matrix(rep(0, nDummySizes * nrow(ism)), ncol = nDummySizes), ism)
            means <- apply(ism, 2, function(x) weighted.mean(0:100, x, na.rm = TRUE))
            ism <- ism[gcPercents, insertSizes]
            heatmap(
                pmin(ism, quantile(ism, 0.975, na.rm = TRUE)), 
                Rowv = NA, Colv = NA,
                scale = "none",
                xlab = "Insert Size (bp)",
                ylab = "Percent GC",
                margins = c(3, 3), 
                labRow = gcPercents, 
                labCol = insertSizes,
                add.expr = points(
                    insertSizes - min(insertSizes) + 1, 
                    means[insertSizes] - min(gcPercents) + 1, 
                    pch = 19,
                    cex = 0.25,
                    col = CONSTANTS$plotlyColors$green
                )
            )
            stopSpinner(session)
        }
    )
    plot
}
insertSizeMatrixPlot("genome")
insertSizeMatrixPlot("spike_in")

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
