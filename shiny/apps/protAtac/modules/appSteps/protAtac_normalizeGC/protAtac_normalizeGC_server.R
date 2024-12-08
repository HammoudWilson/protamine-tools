#----------------------------------------------------------------------
# server components for the protAtac_normalizeGC appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
protAtac_normalizeGCServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'protAtac_normalizeGC'
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
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single") 
samples <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    paBinData(sourceId)$samples
})

#----------------------------------------------------------------------
# samples table
#----------------------------------------------------------------------
sampleTable <- bufferedTableServer(
    "sample",
    id,
    input,
    tableData = reactive( samples()[, .(sample_name, stage)] ),
    selection = 'single',
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
)
sample <- reactive({
    row <- sampleTable$rows_selected()
    req(row)
    samples()[row]
})

#----------------------------------------------------------------------
# appStep outcomes, saved to disk since these are ~one-time analysis steps
# includes GC bias fit, chromosome-level data and junction fits, but not HMM
#----------------------------------------------------------------------
gcBiasFileName <- "gcBiasModels.rds"
gcBiasFile <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    expandSourceFilePath(sourceId, gcBiasFileName)
})
invalidateGcBiasModels <- reactiveVal(1)
getGcBiasModels <- function(gcBiasFile = NULL, sourceId = NULL){
    if(is.null(gcBiasFile)) gcBiasFile <- expandSourceFilePath(sourceId, gcBiasFileName)
    if(file.exists(gcBiasFile)) readRDS(gcBiasFile) else list()
}
gcBiasModels <- reactive({ # gc bias model from negative binomial are calculated synchronously
    invalidateGcBiasModels()
    getGcBiasModels( gcBiasFile = gcBiasFile() )
})
gcBiasModel <- reactive({
    gcBiasModels <- gcBiasModels()
    sample <- sample()
    gcBiasModels[[sample$sample_name]]
})
suppressGcOutliers <- function(nb, gc){
    allowedGC <- range(nb$model$fractionGC)
    pmax(allowedGC[1], pmin(allowedGC[2], gc))
}

#----------------------------------------------------------------------
# interactive GC bias plots, selection cascades to solving negative binomial
#----------------------------------------------------------------------
gcBiasPlotData <- function(){
    sourceId <- sourceId()
    sample <- sample()
    req(sourceId, sample)
    startSpinner(session, message = paste("plotting", sample$sample_name))
    bd <- paBinData(sourceId)
    I <- bd$bins$genome[, excluded == 0 & nAlleles == 2] # 
    nAlleles <- bd$bins$genome[I, nAlleles]
    data.table(
        x = bd$bins$genome[I, pct_gc],
        y = rowSums(bd$binCounts$genome[I, , sample$sample_name], na.rm = TRUE), # rpb = reads per bin
        nAlleles = nAlleles,
        color = ifelse(nAlleles == 2, CONSTANTS$plotlyColors$blue, CONSTANTS$plotlyColors$orange)
    )[sample.int(.N, min(.N, input$nBiasBins))]
}
gcOverplotData <- function(){
    gcBiasModel <- gcBiasModel()
    if(!isTruthy(gcBiasModel)) return(NULL)
    sample <- sample()
    startSpinner(session, message = paste("overplotting", sample$sample_name))
    nb <- gcBiasModel$fit
    gc <- nb$model$fractionGC
    data.table(
        x = gc,
        y = predict(nb, gc, type = 'mu') * 2 # rpba = reads per bin per allele * nAlleles = rpb
    )
}
gcBiasPlot <- interactiveScatterplotServer(
    "gcBiasPlot",
    plotData = reactive({ 
        x <- gcBiasPlotData() 
        stopSpinner(session)
        x
    }),
    accelerate = TRUE,
    # color = reactive({
    #     x <- gcBiasPlotData() 
        
    # }),
    overplot = reactive({
        x <- gcOverplotData()
        stopSpinner(session)
        x
    }),
    overplotMode = "lines",
    overplotColor = CONSTANTS$plotlyColors$red,
    xtitle = "Fraction GC",
    xrange = gcLimits,
    ytitle = "Reads Per Bin (RPB)",
    yrange = function(...) range_pos(..., foldIQR = 5),
    selectable = "lasso"
)
observeEvent(gcBiasPlot$selected(), {
    showUserDialog(
        "Use this GC Bias Fit?", 
        tags$p(
            "If you are happy with your GC bias selection, click OK to run the next, slow actions."
        ), 
        tags$p(
            "If you are not happy, click Cancel and repeat the GC selection."
        ), 
        callback = function(parentInput) {
            removeModal()
            fitGCBiasFromSelected(gcBiasPlot$selected())
        },
        type = 'okCancel', 
        easyClose = FALSE, 
        fade = 250
    )
})
#----------------------------------------------------------------------
gcResidualBiasPlotData <- function(){
    gcBiasModel <- gcBiasModel()
    sourceId <- sourceId()
    sample <- sample()
    req(gcBiasModel, sourceId, sample)
    startSpinner(session, message = paste("plotting", sample$sample_name))
    bd <- paBinData(sourceId)
    I  <- bd$bins$genome[, excluded == 0]
    gc <- bd$bins$genome[I, pct_gc]
    nAlleles <- bd$bins$genome[I, nAlleles]
    rpb_fit <- predict(gcBiasModel$fit, fractionGC = gc, type = "mu") * nAlleles
    data.table(
        x = gc,
        y = (rowSums(bd$binCounts$genome[I, , sample$sample_name], na.rm = TRUE) - rpb_fit) / sqrt(rpb_fit), # rbp z-score under Poisson variance = mean,
        color = ifelse(nAlleles == 2, CONSTANTS$plotlyColors$blue, CONSTANTS$plotlyColors$orange)
    )[sample.int(.N, min(.N, input$nResidualBiasBins))]
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
    ytitle = "Residual RPB Z-Score",
    yrange = function(...) range_both(..., foldIQR = 5),
    selectable = "lasso"
)

#----------------------------------------------------------------------
# fit negative binomial distribution to GC bias
#----------------------------------------------------------------------
fitGCBiasFromSelected <- function(selected){
    req(selected, nrow(selected) > 10)
    sample <- sample()
    startSpinner(session, message = paste("fitting", sample$sample_name))
    selected <- as.data.table(selected)
    gcBiasModels <- gcBiasModels()
    d <- gcBiasPlotData() 
    fit <- new_nbinomCountsGC(
        binCounts  = selected$y,
        fractionGC = selected$x,
        binCN      = d[selected$pointNumber, nAlleles],
        method = 'cubic'
    )
    gcBiasModels[[sample$sample_name]] <- list(
        sample          = sample,
        fit             = fit
    )
    saveRDS(gcBiasModels, file = gcBiasFile())
    stopSpinner(session)
    invalidateGcBiasModels( invalidateGcBiasModels() + 1 ) 
}

#----------------------------------------------------------------------
# external utilities for browser support, etc.
#----------------------------------------------------------------------
getGcBiasModels_externalCall <- function(sourceId){
    gcSourceId <- tryCatch(sourceId(), error = function(e) NULL)
    if(isTruthy(gcSourceId) && gcSourceId == sourceId) gcBiasModels()
    else getGcBiasModels(sourceId = sourceId)
}
getExpectedReadsPerBin <- function(sourceId, sampleName, gc, nAlleles){ # for plotting GC-normalized copy number at bin level
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    gcBiasModel <- gcBiasModels[[sampleName]]
    if(!isTruthy(gcBiasModel)) return(NULL)
    predict(gcBiasModel$fit, suppressGcOutliers(gcBiasModel$fit, gc), type = 'mu') * nAlleles
}
getBinZScore <- function(sourceId, sampleName, gc, nAlleles, binCounts){
    rpb <- getExpectedReadsPerBin(sourceId, sampleName, gc, nAlleles)
    if(!isTruthy(rpb)) return(NULL)
    (binCounts - rpb) / sqrt(rpb) # assuming Poisson variance
}

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    getGcBiasModels_externalCall = getGcBiasModels_externalCall,
    getExpectedReadsPerBin = getExpectedReadsPerBin,
    getBinZScore = getBinZScore,
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
