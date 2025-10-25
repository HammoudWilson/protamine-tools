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
    selection = "multiple"
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
# insert size distribution plots
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
        isdt_this[[series]] * (nInserts_this / nInserts_ref)
    })) %>% as.data.table
    setnames(x, seriesNames)
    x
}
getInsertSize_peakType <- function(n_peak_is_gc_smp, peakType, samples){
    # return a data.table of all insert sizes from 1 to MAX_INSERT_SIZE, one column per requested sample
    peakI <- if(peakType == "non_peak") 0 else 1
    sapply(samples$sample_name, function(sample_name){
        c(
            rep(0, env()$MIN_INSERT_SIZE - 1),
            rowSums(n_peak_is_gc_smp[[peakI + 1]][,,sample_name])
        )
    }, simplify = FALSE, USE.NAMES = TRUE) %>% as.data.table
}
getInsertSizeData <- function(n_peak_is_gc_smp, peakType, samples, aggregate){
    isdt_this <- getInsertSize_peakType(n_peak_is_gc_smp, peakType, samples)
    if(aggregate) isdt_this <- aggregateInsertSizes(isdt_this, samples)
    isdt_this[, lapply(.SD, function(x) x / sum(x, na.rm = TRUE))]
}
insertSizesPlot <- function(peakType, title){
    plot <- staticPlotBoxServer(
        paste("insertSizesPlot", peakType, sep = "_"),
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
            n_peak_is_gc_smp <- paCollate_load_ram(sourceId, "n_peak_is_gc_smp")
            aggregate <- input$aggregateByStage
            isdt <- getInsertSizeData(n_peak_is_gc_smp, peakType, selectedSamples, aggregate)
            seriesNames <- colnames(isdt)
            colors <- if(aggregate) getStageColors(allSamples, selectedSamples) 
                      else if(length(unique(selectedSamples$stage)) == 1) getSampleColorsBySample(selectedSamples)
                      else getSampleColorsByStage(allSamples, selectedSamples)
            plot$initializeFrame(
                xlim = c(env$MIN_INSERT_SIZE, env$MAX_INSERT_SIZE),
                ylim = c(0, max(isdt, na.rm = TRUE)),
                xlab = "Insert Size (bp)",
                ylab = "Frequency",
                xaxs = "i",
                title = title
            )
            abline(v = env$MAPPABILITY_SIZE_LEVELS - 0.5, col = "grey80")
            abline(v = c(meanNucleosomeFragSize, meanDinucleosomeFragSize), col = CONSTANTS$plotlyColors$brown, lty = 3)
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
peakPlot <- insertSizesPlot("peak", "Accessibility Peak Insert Sizes")
nonPeakPlot <- insertSizesPlot("non_peak", "Off-Peak Insert Sizes")

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        peakPlot$settings$replace(bm$outcomes$peakPlotSettings)
        nonPeakPlot$settings$replace(bm$outcomes$nonPeakPlotSettings)
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
        peakPlotSettings = peakPlot$settings$all_,
        nonPeakPlotSettings = nonPeakPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
