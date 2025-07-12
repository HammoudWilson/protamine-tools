#----------------------------------------------------------------------
# server components for the dinucChainsSelectorBox widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
dinucChainsSelectorBoxServer <- function(id, sourceId) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'dinucChainsSelectorBox'
# settings <- activateMdiHeaderLinks( # uncomment as needed
#     session,
#     url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
#     dir = getAppStepDir(module), # for terminal emulator
#     envir = environment(), # for R console
#     baseDirs = getAppStepDir(module), # for code viewer/editor
#     settings = id, # for step-level settings
#     immediate = TRUE # plus any other arguments passed to settingsServer()
# )

#----------------------------------------------------------------------
# data source parsing
#----------------------------------------------------------------------
stages <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTss_ab_initio(sourceId)$stages
})
observeEvent(stages(), {
    stages <- stages()
    req(stages)
    updateSelectInput(
        session, 
        "index_stage", 
        label = "Index Stage",
        choices = c("overlap_group", stages), 
        selected = "overlap_group"
    )
})
applyUserFilter <- function(intervals) {
    stages <- stages()
    req(stages)
    code <- trimws(input$user_filter)
    if(isTruthy(code)){
        for(stage in stages) code <- gsub(stage, paste0(stage, "_rpkm"), code)
        code <- gsub("delta_rpkm", "delta_RPKM", code)
        code <- gsub("max_rpkm", "max_RPKM", code)
        tryCatch(intervals[eval(parse(text = code))], error = function(e) NULL)
    } else {
        intervals
    }
}
unfilteredIntervals <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTss_dinuc_chains(sourceId, input$index_stage, input$include_unpassed)
})
intervals <- reactive({
    unfilteredIntervals() %>%
    applyUserFilter()
})
passedIntervals <- reactive({ # for things only available when RPKM passed
    sourceId <- sourceId()
    req(sourceId)
    paTss_dinuc_chains(sourceId, input$index_stage, FALSE) %>%
    applyUserFilter()
})

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
list(
    stages    = stages,
    unfilteredIntervals = unfilteredIntervals,
    intervals = intervals,
    passedIntervals = passedIntervals,
    applyUserFilter = applyUserFilter,
    input = input
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
