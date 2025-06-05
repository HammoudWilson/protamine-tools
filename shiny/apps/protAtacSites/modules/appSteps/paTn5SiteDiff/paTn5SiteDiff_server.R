#----------------------------------------------------------------------
# server components for the paTn5SiteDiff appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
paTn5SiteDiffServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'paTn5SiteDiff'
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
dinucColors <- list(
    CG = CONSTANTS$plotlyColors$red,
    GC = CONSTANTS$plotlyColors$orange
)

#----------------------------------------------------------------------
# data sources
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "dataSourceTable", 
    selection = "single"
) 

#----------------------------------------------------------------------
# single sample Tn5 plots
#----------------------------------------------------------------------
tn5DiffPlotData <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    sets <- paTn5Sites_tn5SetData(sourceId)
    startSpinner(session, message = "setting MA values")
    I <- sample(sets$n_kmers, 50000)
    list(
        sourceId = sets$sourceId,
        M = sets$f_es[I] - sets$f_rs[I], 
        A = (sets$f_es[I] + sets$f_rs[I]) / 2,
        I = I,
        has_dinuc = sapply(sets$has_dinuc, function(x) x[I], simplify = FALSE, USE.NAMES = TRUE)
    )
})
tn5DiffPlot <- function(dinuc) {
    plot <- staticPlotBoxServer(
        paste("tn5DiffPlot", dinuc, sep = "_"),
        maxHeight = "400px",
        create = function() {
            d <- tn5DiffPlotData()
            dinucColors <- ifelse(
                d$has_dinuc[[dinuc]],
                dinucColors[[dinuc]],
                CONSTANTS$plotlyColors$grey
            ) #  %>% addAlphaToColors(0.25)
            startSpinner(session, message = "rendering diff plot")
            par(mar = c(4.1, 4.1, 2.1, 0.1))
            plot$initializeFrame(
                xlim = c(-8.5,-3.5),
                ylim = c(-1,0.5),
                xlab = "Avg. log10 Site Freq.",
                ylab = "log10 Fold-Change Site Freq. (ES - RS)",
                title = paste("Tn5 MA plot,", dinuc, "dinucleotides")
            )
            abline(h = 0, col = CONSTANTS$plotlyColors$grey)
            plot$addPoints(
                x = d$A,
                y = d$M,
                pch = 16,
                cex = 0.25,
                col = dinucColors
            )
            stopSpinner(session)
        }
    )
}
tn5DiffPlot_CG <- tn5DiffPlot("CG")
tn5DiffPlot_GC <- tn5DiffPlot("GC")

#----------------------------------------------------------------------
# single sample Tn5 plots
#----------------------------------------------------------------------
dinucEnrichmentPlot <- staticPlotBoxServer(
    "dinucEnrichmentPlot",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        d <- paTn5Sites_tn5SetData(sourceId)
        startSpinner(session, message = "rendering violins")
        par(mar = c(4.1, 4.1, 1.1, 0.1))
        dinucEnrichmentPlot$initializeFrame(
            xlim = c(0.5, 6.5),
            ylim = c(-4, 2.5),
            xlab = "",
            ylab = "log10 Enrichment (obs / exp)",
            xaxt = "n"
        )
        axis(1, at = 1:6, labels = c(
            "-", "2+CG", "2+GC",
            "-", "2+CG", "2+GC"
        ))
        mtext("Round", side = 1, line = 2.5, at = 2)
        mtext("Elongating", side = 1, line = 2.5, at = 5)
        abline(h = 0, col = CONSTANTS$plotlyColors$black)
        e_rs_neither <- d$e_rs[d$has_dinuc$neither]
        e_rs_CG      <- d$e_rs[d$has_dinuc$CG]
        e_rs_GC      <- d$e_rs[d$has_dinuc$GC]
        e_es_neither <- d$e_es[d$has_dinuc$neither]
        e_es_CG      <- d$e_es[d$has_dinuc$CG]
        e_es_GC      <- d$e_es[d$has_dinuc$GC]
        vioplot::vioplot(
            e_rs_neither[is.finite(e_rs_neither)],
            e_rs_CG[is.finite(e_rs_CG)],
            e_rs_GC[is.finite(e_rs_GC)],
            e_es_neither[is.finite(e_es_neither)],
            e_es_CG[is.finite(e_es_CG)],
            e_es_GC[is.finite(e_es_GC)],
            ylim = c(-4, 2.5),
            col = c(
                CONSTANTS$plotlyColors$grey,
                dinucColors$CG,
                dinucColors$GC,
                CONSTANTS$plotlyColors$grey,
                dinucColors$CG,
                dinucColors$GC
            ),
            drawRect = TRUE,
            colMed = NA,
            add = TRUE
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
