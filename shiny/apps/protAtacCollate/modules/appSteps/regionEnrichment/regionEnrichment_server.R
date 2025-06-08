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
# BED file upload and display
#----------------------------------------------------------------------
inactivateBedFilesList <- reactiveVal(0)
bedFilesDir <- "uploaded_bed_files"
sourceBedFilesDir <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    expandSourceFilePath(sourceId, bedFilesDir)
})
observeEvent(input$bedFileUpload, {
    sourceId <- tryCatch({
        sourceId()
    }, error = function(e){
        showUserDialog(
            "Data Source Required",
            tags$p("Please select a source data package before uploading BED files."),
            tags$p("BED files are specific to a data package and must be uploaded for each source."),
            callback = function(parentInput) NULL,
            size = "s", 
            type = 'okOnly', 
            easyClose = FALSE, 
            fade = TRUE
        )
        req(FALSE)
    })
    bedFile <- input$bedFileUpload
    bedDir <- sourceBedFilesDir()
    req(bedFile, bedDir)
    bedPath <- file.path(bedDir, bedFile$name)
    file.copy(bedFile$datapath, bedPath, overwrite = TRUE)
    inactivateBedFilesList(inactivateBedFilesList() + 1)
})
bedFileTableData <- reactive({
    inactivateBedFilesList()
    bedDir <- sourceBedFilesDir()
    req(bedDir)
    if(!dir.exists(bedDir)) dir.create(bedDir, recursive = FALSE)
    data.table(
        file = list.files(bedDir, full.names = FALSE)
    )
})
bedFileTable <- bufferedTableServer(
    "bedFileTable",
    id,
    input,  
    bedFileTableData,
    selection = 'single'
)
bedFileName <- reactive({
    row <- bedFileTable$rows_selected()
    req(row)
    bedFileTableData()[row, file]
})
bedData <- reactive({
    row <- bedFileTable$rows_selected()
    bedDir <- sourceBedFilesDir()
    req(row, bedDir)
    bedFile <- file.path(bedDir, bedFileTableData()[row, file])
    bedData <- fread(bedFile, sep = "\t")[, 1:3]
    bed3Cols <- c("chrom", "start0", "end1")
    setnames(bedData, bed3Cols)
    setkeyv(bedData, bed3Cols)
    bedData
})

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
enrichmentPlotData <- reactive({
    sourceId <- sourceId()
    bedData <- bedData()
    req(sourceId, bedData)
    startSpinner(session, message = "getting score data")
    metadata <- paScores_metadata(sourceId)
    config <- list(
        sourceId = sourceId,
        bedFileName = bedFileName(),
        bedData = bedData,
        Scores_Dir = settings$get("Enrichment_Data","Scores_Dir")
    )
    getSampleScores_regions(metadata, config, input$enrichmentScoreType, "score")
})
enrichmentPlot <- function(plot, columns, message){
    nColumns <- length(columns)
    req(nColumns > 0)
    d <- enrichmentPlotData()
    startSpinner(session, message = message)
    scoreType <- scoreTypes$sample[[input$enrichmentScoreType]]
    N_Downsample_Bins <- settings$get("Enrichment_Data","N_Downsample_Bins")
    titleSuffix <- ""
    I <- if(N_Downsample_Bins > 0){
        titleSuffix <- paste0(" (", format(N_Downsample_Bins, big.mark = ",", scientific = FALSE), " total bins)")
        sample(1:nrow(d), N_Downsample_Bins)
    } else TRUE

    par(mar = c(4,4,2,4) + 0.1)
    plot$initializeFrame(
        xlim = c(0.4, nColumns + 0.6),
        ylim = scoreType$valueLim,
        xlab = "", # if(isDelta) paste(label, "Delta") else label,
        ylab = scoreType$enrichmentLabel,
        title = paste(input$enrichmentScoreType, titleSuffix),
        xaxt = "n",
        xaxs = "i"
    )
    plot$addMarginLegend(
        nColumns + 0.7, scoreType$valueLim[2], lty = 1, lwd = 2, 
        legend = c("In", "Out"), bty = "n",
        col = c(CONSTANTS$plotlyColors$green, CONSTANTS$plotlyColors$red), cex = 0.9
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
    colors <- c(
        CONSTANTS$plotlyColors$green,
        CONSTANTS$plotlyColors$red
    )
    colors_ <- addAlphaToColors(colors, 0.5)
    offset <- 0.2
    medians <- list(in_ = numeric(nColumns), out_ = numeric(nColumns))
    for(i in 1:nColumns){
        medians$in_[i]  <- median(d[[columns[i]]][d$hasOverlap == TRUE],  na.rm = TRUE)
        medians$out_[i] <- median(d[[columns[i]]][d$hasOverlap == FALSE], na.rm = TRUE)
        vioplot::vioplot(
            d[I][hasOverlap == TRUE][[columns[i]]], 
            d[I][hasOverlap == FALSE][[columns[i]]], 
            ylim = scoreType$valueLim, 
            names = "in",  
            col = colors_, 
            colMed = colors_, 
            at = i + c(-offset, offset), # placement of pairs on the x-axis
            add = TRUE,
            wex = 0.45, # width expansion factor of the violins
            h = 0.25 # kernel density in # of sd
        )
    }
    lines(
        1:nColumns - offset, 
        medians$in_,
        col = colors[1],
        lwd = 2
    )
    lines(
        1:nColumns + offset, 
        medians$out_,
        col = colors[2],
        lwd = 2
    )
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
