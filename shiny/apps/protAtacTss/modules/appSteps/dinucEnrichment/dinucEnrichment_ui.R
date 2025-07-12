#----------------------------------------------------------------------
# UI components for the dinucEnrichment appStep module
#----------------------------------------------------------------------

# module ui function
dinucEnrichmentUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$dinucEnrichment)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        # settings = TRUE,

        # data source selectors
        fluidRow(
            column(
                width = 8,
                style = "margin: 0; padding: 0;",
                dataSourceTableUI(
                    ns("source"), 
                    "Data Source", 
                    width = 6, 
                    collapsible = FALSE,
                    inFluidRow = FALSE
                ),
                dinucRegionsSelectorBoxUI(
                    ns("dinucRegions"), 
                    width = 6, 
                    includeUmap = FALSE
                ),
                staticPlotBoxUI(
                    ns("overlapPlotBox"), 
                    "Overlap Enrichment",
                    width = 6,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = FALSE,
                    solidHeader = TRUE
                ),
                mdiInteractivePlotBoxUI(
                    ns("scorePlotBox"), 
                    "Score Enrichment",
                    width = 6,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = FALSE,
                    solidHeader = TRUE
                ),
                NULL
            ),
            column(
                width = 4,
                style = "margin: 0; padding: 0;",
                regionsBedTableUI(
                    ns("regionsBedTable"),
                    width = 12
                ),
                clusterProfilePlotBoxUI(
                    ns("clusterProfilePlotBox"), 
                    width = 12
                ),
                NULL
            ),
            NULL
        ),
        regionExpansionUI(
            ns("regionExpansion")
        ),
        NULL
    )
}
