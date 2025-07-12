#----------------------------------------------------------------------
# server components for the dinucRegionsSelectorBox widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
dinucRegionsSelectorBoxServer <- function(id, sourceId) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'dinucRegionsSelectorBox'
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
# data source to region parsing
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
applyUserRegionFilter <- function(regions) {
    stages <- stages()
    req(stages)
    code <- trimws(input$user_filter)
    if(isTruthy(code)){
        for(stage in stages) code <- gsub(stage, paste0(stage, "_rpkm"), code)
        code <- gsub("delta_rpkm", "delta_RPKM", code)
        code <- gsub("max_rpkm", "max_RPKM", code)
        tryCatch(regions[eval(parse(text = code))], error = function(e) NULL)
    } else {
        regions
    }
}
regions <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paTss_dinuc_regions(sourceId, input$index_stage, input$include_unpassed_regions) %>%
    applyUserRegionFilter()
})
passedRegions <- reactive({ # for things only available for passed regions
    sourceId <- sourceId()
    req(sourceId)
    paTss_dinuc_regions(sourceId, input$index_stage, FALSE) %>%
    applyUserRegionFilter()
})

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
list(
    stages  = stages,
    regions = regions,
    passedRegions = passedRegions,
    input = input
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
