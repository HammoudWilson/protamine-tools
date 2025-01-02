#----------------------------------------------------------------------
# server components for the protAtac_vPlots appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
protAtac_vPlotsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'protAtac_vPlots'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    size = "m"
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
samples <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    paTssFragsData(sourceId)$samples
})
tss <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    txnState <- settings$get("V_Plot","Transcription_State")
    tss <- paTssFragsData(sourceId)$tss[[txnState]]
    tss[, tss_I1 := 1:.N]
    tss
})
filteredTssI1 <- reactive({
    tss <- tss()
    req(tss)
    txnState <- settings$get("V_Plot","Transcription_State")

    dstr(tss)
    dstr(txnState)
    dstr(settings$get("V_Plot","Min_RPKM_Active"))
    dstr(settings$get("V_Plot","Min_Fold_Increase"))

    if(txnState == "active"){
        tss[
            gene_start_rpkm   >= settings$get("V_Plot","Min_RPKM_Active") & 
            tss_fold_increase >= settings$get("V_Plot","Min_Fold_Increase"), 
            tss_I1
        ]
    } else {
        maxRpkm <- settings$get("V_Plot","Max_RPKM_Inactive")
        tss[
            upstream_rpkm     <= maxRpkm & 
            gene_start_rpkm   <= maxRpkm & 
            gene_start_rpkm - upstream_rpkm <= maxRpkm / 5, 
            tss_I1
        ]
    }
})

#----------------------------------------------------------------------
# sample selection tables
#----------------------------------------------------------------------
sampleTable <- function(tableId) bufferedTableServer(
    tableId,
    id,
    input,
    tableData = reactive( samples()[, .(sample_name, stage)] ),
    selection = 'multiple',
    options = list(
        paging = FALSE,
        searching = FALSE
    )
)
samples1 <- sampleTable("samples1")
samples2 <- sampleTable("samples2")
getSelectedSamples <- function(bt){
    rows <- bt$rows_selected()
    req(rows)
    samples()[rows]
}
selectedSamples1 <- reactive({ getSelectedSamples(samples1) })
selectedSamples2 <- reactive({ getSelectedSamples(samples2) })

#----------------------------------------------------------------------
# samples V plots
#----------------------------------------------------------------------
binTssFrags <- function(selectedSamplesReactive){
    sourceId <- sourceId()
    req(sourceId)
    selectedSamples <- selectedSamplesReactive()
    req(selectedSamples)
    startSpinner(session, message = "binning frags")
    inc <- settings$get("V_Plot","BP_Per_Bin")
    txnState <- settings$get("V_Plot","Transcription_State")
    lapply(paTssFragsData(sourceId)$tssFrags[[txnState]][selectedSamples$sample_name], function(tf) {
        tf[, ":="(
            x = round((start_to_tss + (end_to_tss - start_to_tss) / 2) / inc, 0) * inc, 
            y = round((end1 - start0) / inc, 0) * inc
        )]
        tf[ , .(tss_i1, x, y)]
    })
}
binnedFrags1 <- reactive({ binTssFrags(selectedSamples1) })
binnedFrags2 <- reactive({ binTssFrags(selectedSamples2) })
#----------------------------------------------------------------------
filterTssFrags <- function(tssFrags){# tssFrags is a list of data.tables, one per sample
    nSamples <- length(tssFrags)
    Max_Distance <- settings$get("V_Plot","Max_X_Distance")
    Max_Size     <- settings$get("V_Plot","Max_Y_Size")

    startSpinner(session, message = "filtering frags")
    filteredTssI1 <- filteredTssI1()
    tssFrags <- lapply(tssFrags, function(tf) {
        tf[
            tss_i1 %in% filteredTssI1 & 
            x >= -Max_Distance & x <= Max_Distance &
            y <= Max_Size,
            .(x, y)
        ]
    })

    startSpinner(session, message = "downsampling frags")
    Max_Plotted_Inserts <- min(
        sapply(tssFrags, nrow), 
        as.integer(settings$get("V_Plot","Max_Plotted_Inserts") / nSamples)
    )
    do.call(rbind, lapply(tssFrags, function(tf) { # each selected sample contributes equally to the plot
        if(nrow(tf) > Max_Plotted_Inserts) tf[sample(.N, Max_Plotted_Inserts)]
        else tf
    }))
}
tssFrags1 <- reactive({ filterTssFrags(binnedFrags1()) })
tssFrags2 <- reactive({ filterTssFrags(binnedFrags2()) })

initializeVPlotFrame <- function(plot, title){
    startSpinner(session, message = "initializing V plot")
    Max_Distance <- settings$get("V_Plot","Max_X_Distance")
    Max_Size     <- settings$get("V_Plot","Max_Y_Size")
    xlim <- c(-Max_Distance, Max_Distance)
    ylim <- c(0, Max_Size)
    par(mar = c(4, 4, 2, 6) + 0.1)
    userTitle <- plot$settings$get("Plot_Frame","Title")
    plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = "Distance to TSS (bp)",
        ylab = "Fragment Size (bp)",
        cex.main = 0.95,
        title = if(!is.null(userTitle) && userTitle != "") userTitle else title
    )
    list(
        x = xlim,
        y = ylim
    )
}
renderVPlot <- function(plot, lim, dt1, dt2 = NULL, h = 210, v = c(-175, 0, 135)){
    startSpinner(session, message = "rendering V plot")
    inc <- settings$get("V_Plot","BP_Per_Bin")
    mdiLevelPlot(
        dt  = dt1,
        dt2 = dt2,
        xlim = lim$x,
        xinc = inc, 
        ylim = lim$y,
        yinc = inc,
        z.fn = length, 
        z.column = "x",
        settings = plot$settings,
        legendTitle = "Count",
        h = h,
        v = v,
        border = NA
    )
}
vPlotServer <- function(plotId, title, tf1Reactive, tf2Reactive = NULL){
    plot <- staticPlotBoxServer(
        plotId,
        title = TRUE,
        settings = mdiLevelPlotSettings,
        size = "m",
        create = function() {
            tf1 <- tf1Reactive()
            req(tf1)
            tf2 <- if(!is.null(tf2Reactive)){
                tf2 <- tf2Reactive()
                req(tf2)
                tf2
                if(nrow(tf1) != nrow(tf2)){
                    stopSpinner(session)
                    message()
                    message("Cannot plot difference: Sample(s) #1 and Sample(s) #2 have different numbers of fragments")
                    message()
                    req(FALSE)
                }
                tf2
            } else NULL
            lim <- initializeVPlotFrame(plot, title)
            renderVPlot(
                plot, 
                lim, 
                dt1 = tf1,
                dt2 = tf2
            )
            stopSpinner(session)
        }
    )
    plot
}
vPlot1     <- vPlotServer("vPlot1",     "Sample(s) #1 V Plot", tssFrags1)
vPlot2     <- vPlotServer("vPlot2",     "Sample(s) #2 V Plot", tssFrags2)
vPlot_diff <- vPlotServer("vPlot_diff", "Sample(s) #2 - Sample(s) #1 V Plot", tssFrags1, tssFrags2)

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
