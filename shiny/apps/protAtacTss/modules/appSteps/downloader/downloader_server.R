#----------------------------------------------------------------------
# server components for the downloader appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
downloaderServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'downloader'
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
observeEvent(input$writeData, {

    # parse share information
    outputDir <- trimws(input$outputDir)
    req(outputDir != "", dir.exists(outputDir))
    sourceId <- sourceId()
    req(sourceId)

    # collect tss ab initio calls
    if(input$scoreType == "dinuc_chains") {
        d <- paTss_ab_initio(sourceId)$intervals

    # collect collate bin scores
    } else {
        startSpinner(session, message = "loading bins")
        bins <- paScores_bins(sourceId)
        startSpinner(session, message = "loading scores")
        scores <- getSampleScores_all(sourceId, input$scoreType)
        startSpinner(session, message = "creating output table")
        d <- merge(
            bins[, .(chrom, start0, end1, included, gc, gc_z, txn)],
            scores,
            by = c("chrom", "start0", "end1"),
            all.x = TRUE,
            sort = FALSE
        )
        rm(bins, scores)
    }

    # report the data structure to web server log
    str(d)

    # write the output table
    startSpinner(session, message = "writing output table")
    if(input$fileFormat == "rds") {
        file <- file.path(outputDir, paste("protAtacTss", input$scoreType, "rds", sep = "."))
        print(file)
        saveRDS(d, file = file)
    } else {
        file <- file.path(outputDir, paste("protAtacTss", input$scoreType, "bed.gz", sep = "."))
        print(file)
        fwrite(
            d, 
            file = file,
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = input$includeHeader
        )
    }
    stopSpinner(session)
})

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
