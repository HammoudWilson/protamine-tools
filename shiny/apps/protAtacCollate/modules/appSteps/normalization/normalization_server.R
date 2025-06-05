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
    selection = "single",
    n_ref_wgt_is_gc_smp = paCollate_load_ram_reactive(sourceId, "n_ref_wgt_is_gc_smp")
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
    selectedSampleName <- selectedSampleName()
    gcBiasModels[[selectedSampleName]]
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
    rpb    <- paCollate_rpb_smp(sourceId, input$normalizeTn5Site, input$normalizeMappability)[, selectedSampleName]
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
    selectedSample <- selectedSample()
    startSpinner(session, message = paste("overplotting", selectedSample$sample_name))
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
# observeEvent(gcBiasPlot$selected(), {
#     showUserDialog(
#         "Use this GC Bias Fit?", 
#         tags$p(
#             "If you are happy with your GC bias selection, click OK to run the next, slow actions."
#         ), 
#         tags$p(
#             "If you are not happy, click Cancel and repeat the GC selection."
#         ), 
#         callback = function(parentInput) {
#             removeModal()
#             fitGCBiasFromSelected(gcBiasPlot$selected())
#         },
#         type = 'okCancel', 
#         easyClose = FALSE, 
#         fade = 250
#     )
# })
#----------------------------------------------------------------------
gcResidualBiasPlotData <- function(){
    gcBiasModel <- gcBiasModel()
    sourceId <- sourceId()
    selectedSampleName <- selectedSampleName()
    req(gcBiasModel, sourceId, selectedSampleName)
    startSpinner(session, message = paste("GC bias", selectedSampleName))
    I <- randomPrimaryBinI()
    gc_bin <- paCollate_gc_bin_smp(sourceId, selectedSampleName) # fractionGC
    rpb    <- paCollate_rpb_smp(sourceId, input$normalizeTn5Site, input$normalizeMappability)[, selectedSampleName]
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
# correlation of Tn5 site weighting vs. GC content
#----------------------------------------------------------------------
tn5vsGCPlotData <- reactive({
    sourceId <- sourceId()
    selectedSampleName <- selectedSampleName()
    req(sourceId, selectedSampleName)
    startSpinner(session, message = paste("plotting", selectedSampleName))
    I <- randomPrimaryBinI()
    gc_bin <- paCollate_gc_bin_smp(sourceId, selectedSampleName) # fractionGC
    n_obs_bin_smp <- paCollate_load_ram(sourceId, "n_obs_bin_smp")[, selectedSampleName]
    n_wgt_bin_smp <- paCollate_load_ram(sourceId, "n_wgt_bin_smp")[, selectedSampleName]
    data.table(
        x = gc_bin[I], # same as fractionGC
        y = log10(n_wgt_bin_smp[I] / n_obs_bin_smp[I]) # Tn5 site weighting
    )
})
tn5vsGCPlot <- staticPlotBoxServer(
    "tn5vsGCPlot",
    maxHeight = "400px",
    create = function() {
        d <- tn5vsGCPlotData()
        fit <- lm(y ~ x, data = d)
        par(mar = c(4, 4, 0, 0) + 0.1)
        tn5vsGCPlot$initializeFrame(
            xlim = gcLimits,
            ylim = c(-3,3),
            xlab = "Fraction GC",
            ylab = "log10(n_wgt / n_obs)"
        )
        abline(h = -5:5, col = CONSTANTS$plotlyColors$grey)
        tn5vsGCPlot$addPoints(
            x = d$x,
            y = d$y,
            pch = 16,
            cex = 0.25,
            col = CONSTANTS$plotlyColors$blue
        )
        abline(
            fit,
            col = CONSTANTS$plotlyColors$red,
            lwd = 1.5
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# examine Tn5 weight distribution per insert
#----------------------------------------------------------------------
wgtDistPlot <- staticPlotBoxServer(
    "wgtDistPlot",
    maxHeight = "400px",
        create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        req(sourceId, selectedSampleName)
        d <- paCollate_n_ins_wgt_smp(sourceId, selectedSampleName)
        d$log10_wgt <- d$log10_wgt / 10 # since value comes to us as int(round(log10(x) * 10))
        q <- quantile(rep(10 ** d$log10_wgt, d$n_ins), c(0.01, 0.5, 0.99)) %>% log10()
        par(mar = c(4, 4, 0, 0) + 0.1)
        wgtDistPlot$initializeFrame(
            xlim = c(-5, 5),
            ylim = c(0, max(d$n_ins) * 1.05),
            xlab = "log10 Tn5 Weight",
            ylab = "# of Inserts"
        )
        abline(v = 0, col = CONSTANTS$plotlyColors$black, lty = 1)
        wgtDistPlot$addPoints(
            x = d$log10_wgt,
            y = d$n_ins,
            pch = 16,
            col = CONSTANTS$plotlyColors$blue
        )
        abline(v = q, col = CONSTANTS$plotlyColors$red, lty = 1)
        stopSpinner(session)
    }
)

# #----------------------------------------------------------------------
# # fit negative binomial distribution to GC bias
# #----------------------------------------------------------------------
# fitGCBiasFromSelected <- function(selected){
#     req(selected, nrow(selected) > 10)
#     selectedSample <- selectedSample()
#     startSpinner(session, message = paste("fitting", selectedSample$sample_name))
#     gcBiasModels <- gcBiasModels()

#     sourceId <- sourceId()
#     req(sourceId)
#     bins <- paCollate_bins(sourceId)
#     gc_bin <- paCollate_gc_bin_smp(sourceId, selectedSample$sample_name) # fractionGC
#     rpb    <- paCollate_rpb_smp(sourceId, input$normalizeTn5Site, input$normalizeMappability)[, selectedSample$sample_name]
#     genome <- paCollate_env(sourceId)$PRIMARY_GENOME

#     I <- getIncludedAutosomeBins(bins, gc_bin, genome)

#     dt <- data.table(
#         gc_bin = gc_bin,
#         rpb    = rpb
#     )[I, 
#         .(
#             q025 = quantile(rpb, 0.025),
#             q975 = quantile(rpb, 0.975),
#             gc_bin,
#             rpb
#         ), 
#         keyby = .(x = as.integer(round(gc_bin * 20)))
#     ][rpb >= q025 & rpb <= q975]
#     J <- sample(1:nrow(dt), 200000)
#     fit <- new_nbinomCountsGC(
#         binCounts  = dt$rpb[J],
#         fractionGC = dt$gc_bin[J],
#         binCN      = 2,
#         method = 'cubic'
#     )
#     dt <- dt[zScore(fit, dt$rpb, dt$gc_bin, 2) %between% c(-3, 3)]
#     J <- sample(1:nrow(dt), 200000)
#     fit <- new_nbinomCountsGC(
#         binCounts  = dt$rpb[J],
#         fractionGC = dt$gc_bin[J],
#         binCN      = 2,
#         method = 'cubic'
#     )

#     # selected <- as.data.table(selected)
#     # d <- gcBiasPlotData() 
#     # fit <- new_nbinomCountsGC(
#     #     binCounts  = selected$y,
#     #     fractionGC = selected$x,
#     #     binCN      = 2, # d[selected$pointNumber, nAlleles],
#     #     method = 'cubic'
#     # )

#     gcBiasModels[[selectedSample$sample_name]] <- list(
#         sample = selectedSample,
#         fit    = fit
#     )
#     saveRDS(gcBiasModels, file = gcBiasFile())
#     stopSpinner(session)
#     invalidateGcBiasModels( invalidateGcBiasModels() + 1 ) 
# }

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
