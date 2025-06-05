#----------------------------------------------------------------------
# server components for the paSampleCounts appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
paSampleCountsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'paSampleCounts'
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

#----------------------------------------------------------------------
# sample insert counts table
#----------------------------------------------------------------------
paSampleCountsTableData <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    merge(
        paTn5Sites_samples(sourceId),
        fread(getSourceFilePath(sourceId, "n-inserts-by-stage")),
        by = "sample_name",
        all.x = TRUE,
        sort = FALSE
    )[, .( 
        order = staging_order, 
        batch, 
        prefix = filename_prefix, 
        sample_name,
        staging, 
        stage, 
        proper_pairs = stage_1_proper_pairs, # 1
        included = stage_2_included, # 2
        dedup = stage_3_dedup, # 3
        mappable = stage_4_mappable, # 4
        frac_included = frac_included_2_1, 
        frac_unique = frac_unique_3_2, 
        frac_mappable = frac_mappable_4_3
    )]
})
bufferedTableServer(
    "paSampleCountsTable",
    id,
    input,  
    paSampleCountsTableData,
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
