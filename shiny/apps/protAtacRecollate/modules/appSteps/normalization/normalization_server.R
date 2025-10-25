#----------------------------------------------------------------------
# server components for the normalization appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
normalizationServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'normalization'
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
spermatidSamplesTable <- spermatidSamplesTableServer(
    "spermatidSamplesTable", 
    samples = paCollate_samples_reactive(sourceId),
    selection = "single"
)
selectedSample <- reactive({
    spermatidSamplesTable$selectedSamples()
})
selectedSampleName <- reactive({
    selectedSample()$sample_name
})
env <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paCollate_env(sourceId)
})

#----------------------------------------------------------------------
# appStep outcomes, saved to disk since these are ~one-time analysis steps
# includes GC bias fit, chromosome-level data and junction fits, but not HMM
#----------------------------------------------------------------------
gcBiasModel <- reactive({
    sourceId <- sourceId()
    selectedSampleName <- selectedSampleName()
    req(sourceId, selectedSampleName)
    gcBiasModels <- paCollate_gc_bias_models(sourceId)
    grcz_type <- "gcrz_obs"
    list(
        sample_name = selectedSampleName,
        fit = gcBiasModels[[grcz_type]][[selectedSampleName]]
    )
})

#----------------------------------------------------------------------
# interactive GC bias plots, selection cascades to solving negative binomial
#----------------------------------------------------------------------
randomPrimaryBinI <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "setting random bins")
    bins <- paCollate_bins(sourceId)
    gc_bin <- rowMeans(paCollate_load_ram(sourceId, "gc_bin_smp"))
    genome <- paCollate_env(sourceId)$PRIMARY_GENOME
    I <- getIncludedAutosomeBins(bins, gc_bin, genome)
    stopSpinner(session)
    sample(which(I), settings$get("Dot_Plots","N_Plotted_Bins"))
})
gcBiasPlotData <- reactive({
    sourceId <- sourceId()
    selectedSampleName <- selectedSampleName()
    req(sourceId, selectedSampleName)
    startSpinner(session, message = paste("plotting", selectedSampleName))
    I <- randomPrimaryBinI()
    gc_bin <- paCollate_gc_bin_smp(sourceId, selectedSampleName) # fractionGC
    rpb    <- paCollate_rpb_smp(sourceId, FALSE)[, selectedSampleName]
    data.table(
        x = gc_bin[I], # same as fractionGC
        y = rpb[I],    # rpb = reads per bin
        nAlleles = 2,  # same as binCN
        color = CONSTANTS$plotlyColors$blue
    )
})
gcOverplotData <- reactive({
    gcBiasModel <- gcBiasModel()
    if(!isTruthy(gcBiasModel)) return(NULL)
    startSpinner(session, message = paste("overplotting", gcBiasModel$sample_name))
    nb <- gcBiasModel$fit
    gc <- nb$model$fractionGC
    data.table(
        x = gc,
        y = predict(nb, gc, type = 'mu') * 2 # rpba = reads per bin per allele * nAlleles = rpb
    )
})
gcBiasPlot <- interactiveScatterplotServer(
    "gcBiasPlot",
    plotData = reactive({ 
        x <- gcBiasPlotData() 
        stopSpinner(session)
        x
    }),
    accelerate = TRUE,
    overplot = reactive({
        x <- gcOverplotData()
        stopSpinner(session)
        x
    }),
    overplotMode = "lines",
    overplotColor = CONSTANTS$plotlyColors$red,
    xtitle = "Fraction GC",
    xrange = gcLimits,
    ytitle = "Reads Per Bin",
    yrange = function(...) range_pos(..., foldIQR = 5),
    selectable = FALSE # "lasso"
)
#----------------------------------------------------------------------
gcResidualBiasPlotData <- function(){
    gcBiasModel <- gcBiasModel()
    sourceId <- sourceId()
    req(gcBiasModel, sourceId)
    startSpinner(session, message = paste("GC bias", gcBiasModel$sample_name))
    I <- randomPrimaryBinI()
    gc_bin <- paCollate_gc_bin_smp(sourceId, gcBiasModel$sample_name) # fractionGC
    rpb    <- paCollate_rpb_smp(sourceId, FALSE)[, gcBiasModel$sample_name]
    binCN <- 2 # bd$bins$genome[I, nAlleles]
    data.table(
        x = gc_bin[I],
        y = zScore(gcBiasModel$fit, rpb[I], gc_bin[I], binCN),
        color = CONSTANTS$plotlyColors$blue # ifelse(binCN == 2, CONSTANTS$plotlyColors$blue, CONSTANTS$plotlyColors$orange)
    )
}
gcResidualBiasPlot <- interactiveScatterplotServer(
    "gcResidualBiasPlot",
    plotData = reactive({ 
        x <- gcResidualBiasPlotData()
        stopSpinner(session)
        x
    }),
    accelerate = TRUE,
    color = CONSTANTS$plotlyColors$blue,
    xtitle = "Fraction GC",
    xrange = gcLimits,
    ytitle = "Residual Reads Per Bin Z Score",
    yrange = function(...) range_both(..., foldIQR = 5),
    selectable = "lasso"
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
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
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
