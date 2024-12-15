#----------------------------------------------------------------------
# server components for the protAtac_txnCorrelation appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
protAtac_txnCorrelationServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'protAtac_txnCorrelation'
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
sampleToPlot <- reactive({
    row <- sampleTable$rows_selected()
    req(row)
    samples()[row]
})

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
txnCorrelationPlot <- staticPlotBoxServer(
    "txnCorrelationPlot",
    maxHeight = "400px",
    points   = TRUE,
    legend  = FALSE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        sampleToPlot <- sampleToPlot()
        req(sampleToPlot)
        xScore <- getGenomeScores(sourceId, "txn")[[1]]$score
        yScore <- getSampleScores(sourceId, input$scoreType, sampleToPlot)[[1]]$score
        I <- sample(1:length(xScore), 20000)
        xScore <- xScore[I]
        yScore <- yScore[I]
        fit <- lm(yScore ~ xScore)
        startSpinner(session, message = "plotting")
        txnCorrelationPlot$initializeFrame(
            xlim = c(-2, 2),
            ylim = range(yScore, na.rm = TRUE),
            xlab = paste(scoreTypes$genome$txn$label, scoreTypes$genome$txn$unit),
            ylab = paste(scoreTypes$sample[[input$scoreType]]$label, scoreTypes$sample[[input$scoreType]]$unit)
        )
        txnCorrelationPlot$addPoints(
            x   = xScore,
            y   = yScore,
            col = rgb(0, 0, 0, 0.25),
            pch = 16,
            cex = 0.25
        )
        txnCorrelationPlot$addLine(
            x = c(-2, 2),
            y = c(-2, 2) * coef(fit)[2] + coef(fit)[1],
            col = paColors$BLUE,
            lwd = 1.5
        )
        stopSpinner(session)
    }
)
#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        # genomePlot$settings$replace(bm$outcomes$genomePlotSettings)
        # spikeInPlot$settings$replace(bm$outcomes$spikeInPlotSettings)
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
        # genomePlotSettings  = genomePlot$settings$all_,
        # spikeInPlotSettings = spikeInPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
