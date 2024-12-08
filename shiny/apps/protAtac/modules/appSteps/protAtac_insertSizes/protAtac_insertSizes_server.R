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
    settings = id, # for step-level settings
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
            req(sourceId, samples)
            isd <- paInsertSizes(sourceId)
            aggregate <- settings$get("Insert_Sizes","Aggregate_Samples_By_Stage")
            normalize <- settings$get("Insert_Sizes","Normalize_To_Spike_In")
            if(refType == "spike_in" && normalize) req(FALSE)
            binSize <- isd$bin_size
            isd <- getInsertSizeData(isd$insertSizes, refType, samples, aggregate, normalize)
            seriesNames <- colnames(isd)
            colors <- stageColors[if(aggregate) seriesNames else samples$stage]
            # colors <- 1:length(seriesNames)
            names(colors) <- seriesNames
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
            abline(v = 146, col = CONSTANTS$plotlyColors$grey) # nucleosome size
            plot$addLegend(
                legend = seriesNames,
                col = colors,
                cex = 0.85
            )
        }
    )
    plot
}
genomePlot  <- insertSizesPlot("genome")
spikeInPlot <- insertSizesPlot("spike_in")

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
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
    settings = settings$all_,
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
