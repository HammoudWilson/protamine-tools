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
    samples = paCollate_samples_reactive(sourceId),
    selection = "multiple",
    n_ref_wgt_is_gc_smp = paCollate_load_ram_reactive(sourceId, "n_ref_wgt_is_gc_smp")
)
selectedSamples <- reactive({
    spermatidSamplesTable$selectedSamples()
})
selectedSampleNames <- reactive({
    selectedSamples()$sample_name
})
allSamples <- reactive({
    spermatidSamplesTable$allSamples()
})
env <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paCollate_env(sourceId)
})

#----------------------------------------------------------------------
# insert size vs. GC matrix heatmap plots
#----------------------------------------------------------------------
insertSizeMatrixPlot <- function(refType){
    plot <- staticPlotBoxServer(
        paste("insertSizeMatrix", refType, sep = "_"),
        maxHeight = "400px",
        create = function() {
            sourceId <- sourceId()
            selectedSampleNames <- selectedSampleNames()
            env <- env()
            req(sourceId, selectedSampleNames, env)
            n_is_gc_smp <- paCollate_load_ram(sourceId, "n_ref_wgt_is_gc_smp")[[refType]]$observed
            n_is_gc <- apply(n_is_gc_smp[,,selectedSampleNames], c(1, 2), sum)
            n_is_gc <- if(input$scaleHeatmapByInsertSize){
                apply(n_is_gc, 1, function(x) x / sum(x, na.rm = TRUE))
            } else {
                t(n_is_gc)
            }
            gcPercents <- 25:75 # hardcoded here just for dev display purposes
            minInsertSize <- if(input$nucleosomesOnly) 120 else env$MIN_INSERT_SIZE
            insertSizes <- minInsertSize:env$MAX_INSERT_SIZE
            nDummySizes <- env$MIN_INSERT_SIZE - 1
            n_is_gc <- cbind(matrix(rep(0, nDummySizes * nrow(n_is_gc)), ncol = nDummySizes), n_is_gc)
            means <- apply(n_is_gc, 2, function(x) weighted.mean(0:100, x, na.rm = TRUE))
            n_is_gc <- n_is_gc[gcPercents, insertSizes]
            heatmap(
                pmin(n_is_gc, quantile(n_is_gc, 0.975, na.rm = TRUE)), 
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
insertSizeMatrixPlot("primary")
insertSizeMatrixPlot("spike_in")

# #----------------------------------------------------------------------
# # insert size distribution plots
# #----------------------------------------------------------------------
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
getInsertSize_refType <- function(n_ref_wgt_is_gc_smp, refType, samples){
    # return a data.table of all insert sizes from 1 to MAX_INSERT_SIZE, one column per requested sample
    sapply(samples$sample_name, function(sample_name){
        c(
            rep(0, env()$MIN_INSERT_SIZE - 1),
            rowSums(n_ref_wgt_is_gc_smp[[refType]]$observed[,,sample_name])
        )
    }, simplify = FALSE, USE.NAMES = TRUE) %>% as.data.table
}
getInsertSizeData <- function(n_ref_wgt_is_gc_smp, refType, samples, aggregate, normalize){
    isdt_this <- getInsertSize_refType(n_ref_wgt_is_gc_smp, refType, samples)
    if(aggregate) isdt_this <- aggregateInsertSizes(isdt_this, samples)
    if(normalize){
        isdt_ref <-getInsertSize_refType(n_ref_wgt_is_gc_smp, "spike_in", samples)
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
            env <- env()
            req(sourceId, selectedSamples, env)
            allSamples <- allSamples()
            n_ref_wgt_is_gc_smp <- paCollate_load_ram(sourceId, "n_ref_wgt_is_gc_smp")
            aggregate <- input$aggregateByStage
            normalize <- input$normalizeToSpikeIn
            if(refType == "spike_in" && normalize) req(FALSE)
            isdt <- getInsertSizeData(n_ref_wgt_is_gc_smp, refType, selectedSamples, aggregate, normalize)
            seriesNames <- colnames(isdt)
            colors <- if(aggregate) getStageColors(allSamples, selectedSamples) 
                      else if(length(unique(selectedSamples$stage)) == 1) getSampleColorsBySample(selectedSamples)
                      else getSampleColorsByStage(allSamples, selectedSamples)
            plot$initializeFrame(
                xlim = c(env$MIN_INSERT_SIZE, env$MAX_INSERT_SIZE),
                ylim = c(0, max(isdt, na.rm = TRUE)),
                xlab = "Insert Size (bp)",
                ylab = "Frequency",
                xaxs = "i"
            )
            abline(v = env$MAPPABILITY_SIZE_LEVELS - 0.5, col = "grey80")
            for(series in seriesNames){
                plot$addLines( # addLines follows the same pattern, etc.
                    x = env$MIN_INSERT_SIZE:env$MAX_INSERT_SIZE,
                    y = isdt[[series]][env$MIN_INSERT_SIZE:env$MAX_INSERT_SIZE],
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
primaryPlot <- insertSizesPlot("primary")
spikeInPlot <- insertSizesPlot("spike_in")

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        primaryPlot$settings$replace(bm$outcomes$primaryPlotSettings)
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
        primaryPlotSettings = primaryPlot$settings$all_,
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
