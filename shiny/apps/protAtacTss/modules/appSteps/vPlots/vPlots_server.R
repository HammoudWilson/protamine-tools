#----------------------------------------------------------------------
# server components for the vPlots appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
vPlotsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'vPlots'
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
# empirically determined limits
#----------------------------------------------------------------------
meanNucDist_proximal   <- -175
meanNucDist_distal     <-  133
DIAG_X1 <- -170 
DIAG_Y1 <- meanNucleosomeFragSize
DIAG_X2 <- 0
DIAG_Y2 <- 125
interNucDist <- -310
minAllowedSize_nuc <- 150
maxAllowedSize <- 450
maxAllowedDist_nuc <- 350
plotDistLim <- c(-maxAllowedSize, maxAllowedDist_nuc)
plotSizeLim <- c(0, maxAllowedSize)
plotNucV <- c(interNucDist, DIAG_X1, 0, meanNucDist_distal)
plotLim <- list(
    TSS = list(
        x = plotDistLim + meanNucDist_distal,
        y = plotSizeLim
    ),
    nucleosome = list(
        x = plotDistLim,
        y = plotSizeLim
    )
)
plotV <- list(
    TSS = c(DIAG_X1, 0, meanNucDist_distal),
    nucleosome = c(interNucDist, DIAG_X1, 0, meanNucDist_distal)
)
nucDistLim <- c(DIAG_X1, maxAllowedDist_nuc)
nucSizeLim <- c(minAllowedSize_nuc, maxAllowedSize)
nucLim <- list(
    TSS = list(
        x = nucDistLim + meanNucDist_distal,
        y = nucSizeLim
    ),
    nucleosome = list(
        x = nucDistLim,
        y = nucSizeLim
    )
)
is_above_line_ <- function(x, y) {
    slope <- (DIAG_Y2 - DIAG_Y1) / (DIAG_X2 - DIAG_X1)
    intercept <- DIAG_Y1 - slope * (DIAG_X1 + meanNucDist_distal)
    y > slope * x + intercept
}

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
samples <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    paTss_frags(sourceId)$samples
})
tss <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    txnState <- settings$get("V_Plot","Transcription_State")
    paTss_frags(sourceId)$tss[[txnState]]
})
filteredTssI1 <- reactive({
    tss <- tss()
    req(tss)
    txnState <- settings$get("V_Plot","Transcription_State")
    plotCenter <- settings$get("V_Plot","Plot_Center")
    if(txnState == "active"){
        tss[
            gene_start_rpkm   >= settings$get("V_Plot","Min_RPKM_Active") & 
            tss_fold_increase >= settings$get("V_Plot","Min_Fold_Increase") & 
            if(plotCenter == "TSS") TRUE else !is.na(nucleosome_distance), 
            tss_i1
        ]
    } else {
        req(plotCenter == "TSS")
        maxRpkm <- settings$get("V_Plot","Max_RPKM_Inactive")
        tss[
            upstream_rpkm     <= maxRpkm & 
            gene_start_rpkm   <= maxRpkm & 
            gene_start_rpkm - upstream_rpkm <= maxRpkm / 5, 
            tss_i1
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
    lapply(paTss_frags(sourceId)$tssFrags[[txnState]][selectedSamples$sample_name], function(tf) {
        tf[, ":="(
            x = round((start_to_tss + (end_to_tss - start_to_tss) / 2) / inc, 0) * inc, # distance before applying nucleosome distance correction
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
    Plot_Center  <- settings$get("V_Plot","Plot_Center")
    Nucleosome_Only <- settings$get("V_Plot","Nucleosome_Only")

    startSpinner(session, message = "filtering frags")
    filteredTssI1 <- filteredTssI1()
    tss <- tss()
    filterLim <- if(Nucleosome_Only) nucLim$TSS else plotLim$TSS # use TSS limits since haven't applied nucleosome distance correction yet
    tssFrags <- lapply(tssFrags, function(tf) {
        tf[
            tss_i1 %in% filteredTssI1 & 
            x %between% filterLim$x &
            y %between% filterLim$y &
            (!Nucleosome_Only | is_above_line_(x, y)),
            {
                if(Plot_Center == "nucleosome") {
                    tss_i1_ <- tss_i1
                    x <- x - tss[tss_i1_, nucleosome_distance]
                }
                .(x, y)
            }
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
    Plot_Center  <- settings$get("V_Plot","Plot_Center")
    par(mar = c(4, 4, 2, 6) + 0.1)
    userTitle <- plot$settings$get("Plot_Frame","Title")
    plot$initializeFrame(
        xlim = plotLim[[Plot_Center]]$x,
        ylim = plotLim[[Plot_Center]]$y,
        xlab = paste("Distance to", Plot_Center, "(bp)"),
        ylab = "Fragment Size (bp)",
        cex.main = 0.95,
        title = if(!is.null(userTitle) && userTitle != "") userTitle else title
    )
    plotLim[[Plot_Center]]
}
renderVPlot <- function(plot, lim, dt1, dt2 = NULL){
    startSpinner(session, message = "rendering V plot")
    inc <- settings$get("V_Plot","BP_Per_Bin")
    Plot_Center  <- settings$get("V_Plot","Plot_Center")
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
        h = c(meanNucleosomeFragSize, meanDinucleosomeFragSize),
        v = plotV[[Plot_Center]],
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
                if(nrow(tf1) > nrow(tf2)){
                    tf1 <- tf1[sample(.N, nrow(tf2))]
                } else if(nrow(tf2) > nrow(tf1)){
                    tf2 <- tf2[sample(.N, nrow(tf1))]
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
