#----------------------------------------------------------------------
# server components for the regionEnrichment appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
regionEnrichmentServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'regionEnrichment'
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
    selection = "multiple",
    n_ref_wgt_is_gc_smp = paCollate_load_ram_reactive(sourceId, "n_ref_wgt_is_gc_smp")
)
selectedSamples <- reactive({
    spermatidSamplesTable$selectedSamples()
})
selectedSampleNames <- reactive({
    selectedSamples()$sample_name
})
selectedStages <- reactive({
    unique(selectedSamples()$stage)
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
# BED file selection
#----------------------------------------------------------------------
regionsBedTable <- regionsBedTableServer(
    "regionsBedTable",
    sourceId
)

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
enrichmentPlotData <- reactive({
    sourceId <- sourceId()
    bedData <- regionsBedTable$data()
    req(sourceId, bedData)
    startSpinner(session, message = "getting score data")
    metadata <- paScores_metadata(sourceId)
    config <- list(
        sourceId = sourceId,
        bedFileName = regionsBedTable$name(),
        bedData = bedData,
        Scores_Dir = settings$get("Enrichment_Data","Scores_Dir")
    )
    dstr(config)
    getSampleScores_regions(metadata, config, input$enrichmentScoreType, "score")
})
enrichmentPlot <- function(plot, columns, message){
    nColumns <- length(columns)
    req(nColumns > 0)
    d <- enrichmentPlotData()
    startSpinner(session, message = message)
    scoreType <- scoreTypes$sample[[input$enrichmentScoreType]]
    N_Downsample_Bins <- settings$get("Enrichment_Data","N_Downsample_Bins")
    titleSuffix <- if(N_Downsample_Bins > 0){
        paste0(" (", format(N_Downsample_Bins, big.mark = ",", scientific = FALSE), " downsample bins)")
    } else ""

    par(mar = c(4,4,2,4) + 0.1)
    ylim <- as.numeric(strsplit(trimws(settings$get("Enrichment_Data","Y_Axis_Limits")), ",")[[1]])
    if(is.na(ylim[2])) ylim <- c(-abs(ylim[1]), abs(ylim[1]))
    plot$initializeFrame(
        xlim = c(0.4, nColumns + 0.6),
        ylim = ylim,
        xlab = "",
        ylab = scoreType$enrichmentLabel,
        title = paste(input$enrichmentScoreType, titleSuffix),
        xaxt = "n",
        xaxs = "i"
    )
    colors <- c(
        CONSTANTS$plotlyColors$grey,
        CONSTANTS$plotlyColors$green,
        CONSTANTS$plotlyColors$red
    )
    plot$addMarginLegend(
        nColumns + 0.7, ylim[2], lty = 1, lwd = 2, 
        legend = c("All", "In", "Out"), bty = "n",
        col = colors, cex = 0.9
    )
    axis(1, at = 1:nColumns, labels = FALSE)
    text(
        x = 1:nColumns, 
        y = par("usr")[3] - 0.1 * diff(par("usr")[3:4]), 
        labels = columns, 
        srt = 35, 
        adj = 1, 
        xpd = TRUE, 
        cex = 0.9
    )
    abline(h = 0)

    colors_ <- addAlphaToColors(colors, 0.5)
    xOffset <- 0.3
    xOffsets <- c(-xOffset, 0, xOffset) # offsets for violin pairs
    medians <- list(
        all_ = numeric(nColumns),
        in_  = numeric(nColumns), 
        out_ = numeric(nColumns)
    )
    downsampleBins <- function(d) {
        n <- length(d)
        if(N_Downsample_Bins > 0 && n > N_Downsample_Bins) d[sample(1:n, N_Downsample_Bins)] else d
    }
    for(columnI in 1:nColumns){
        d_all <- d[[columns[columnI]]]
        d_in  <- d_all[d$hasOverlap == TRUE]
        d_out <- d_all[d$hasOverlap == FALSE]
        medians$all_[columnI] <- median(d_all, na.rm = TRUE)
        medians$in_[columnI]  <- median(d_in,  na.rm = TRUE)
        medians$out_[columnI] <- median(d_out, na.rm = TRUE)
        vioplot::vioplot(
            downsampleBins(d_all), 
            downsampleBins(d_in), 
            downsampleBins(d_out), 
            ylim = ylim, 
            names = "", # not used
            col = colors_, 
            colMed = colors_, 
            at = columnI + xOffsets, # placement of pairs on the x-axis
            add = TRUE,
            wex = 0.33, # width expansion factor of the violins
            h = 0.25 # kernel density in # of sd
        )
    }
    for(series in 1:3){
        lines(
            1:nColumns + xOffsets[series], 
            medians[[series]], 
            col = colors[series],
            lwd = 2
        )
    }
    stopSpinner(session)
}
bySamplePlot <- staticPlotBoxServer(
    "bySamplePlot",
    maxHeight = "400px",
    create = function() {
        columns <- selectedSampleNames()
        req(columns)
        enrichmentPlot(bySamplePlot, columns, message = "rendering sample plot")
    }
)
byStagePlot <- staticPlotBoxServer(
    "byStagePlot",
    maxHeight = "400px",
    create = function() {
        columns <- selectedStages()
        req(columns)
        enrichmentPlot(byStagePlot, columns, message = "rendering stage plot")
    }
)
byStageTypePlot <- staticPlotBoxServer(
    "byStageTypePlot",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        metadata <- paScores_metadata(sourceId)
        stopSpinner(session)
        columns <- names(metadata$stageTypes)
        req(columns)
        enrichmentPlot(byStageTypePlot, columns, message = "rendering stage type plot")
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
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
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
