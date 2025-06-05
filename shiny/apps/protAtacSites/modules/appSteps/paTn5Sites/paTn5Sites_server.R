#----------------------------------------------------------------------
# server components for the paTn5Sites appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
paTn5SitesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'paTn5Sites'
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
spermatidSamplesTable <- spermatidSamplesTableServer(
    "spermatidSamplesTable", 
    samples = paTn5Sites_samples_reactive(sourceId),
    selection = "single"
)
selectedSample <- reactive({
    spermatidSamplesTable$selectedSamples()
})
selectedSampleName <- reactive({
    selectedSample()$sample_name
})

#----------------------------------------------------------------------
# single sample Tn5 plots
#----------------------------------------------------------------------
tn5CdfPlot <- staticPlotBoxServer(
    "tn5CdfPlot",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        req(selectedSampleName)
        startSpinner(session, message = "loading site data")
        ordered_sites <- paTn5Sites_table(sourceId, by = input$tn5Site_rankBy)
        f_prm_obs_smp <- paTn5Sites_f_prm_obs_smp(sourceId, selectedSampleName)[ordered_sites$order]
        cdf_smp <- cumsum(f_prm_obs_smp)
        startSpinner(session, message = "rendering Tn5 cdf plot")
        nSites <- nrow(ordered_sites$dt)
        tn5CdfPlot$initializeFrame(
            xlim = c(1, nSites),
            ylim = c(0, 1),
            xlab = "Tn5 Kmer Rank",
            ylab = "Cumulative Fraction of Tn5 Sites"
        )
        i <- sample(1:nSites, 20000)
        # add this sample's CDF
        tn5CdfPlot$addPoints(
            x = i,
            y = cdf_smp[i],
            pch = 16,
            cex = 0.75,
            col = CONSTANTS$plotlyColors$blue
        )
        # overplot with the observed CDF aggregated over all samples
        tn5CdfPlot$addPoints(
            x = i,
            y = ordered_sites$dt$cdf_obs[i],
            pch = 16,
            cex = 0.35,
            col = CONSTANTS$plotlyColors$grey
        )
        # add the expected CDF from the genome
        tn5CdfPlot$addPoints(
            x = i,
            y = ordered_sites$dt$cdf_exp[i],
            pch = 16,
            cex = 0.75,
            col = CONSTANTS$plotlyColors$green
        )
        # add a reference line for the uniform distribution of site frequencies
        tn5CdfPlot$addLines(
            x = c(1, nSites),
            y = c(0, 1),
            col = CONSTANTS$plotlyColors$grey
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# Tn5 sites table
#----------------------------------------------------------------------
tn5SiteTableData <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading sites table")
    dt <- paTn5Sites_table(sourceId, by = input$tn5Site_rankBy)$dt
    stopSpinner(session)
    req(dt)
    dt[c(1:100, (.N-100):.N), .(
        rank = rank,
        kmer = kmer,
        site = paTn5Sites_expand(kmer),
        log2_e      = e %>% log2  %>% round(2),
        log10_f_obs = f_obs %>% log10 %>% round(2),
        log10_f_exp = f_exp %>% log10 %>% round(2)
    )]
})
bufferedTableServer(
    "tn5SiteTable",
    id,
    input,  
    tn5SiteTableData,
    selection = 'none'
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
