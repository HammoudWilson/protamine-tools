#----------------------------------------------------------------------
# UI components for the dinucRegions appStep module
#----------------------------------------------------------------------

# module ui function
dinucRegionsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$dinucRegions)

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
            dataSourceTableUI(
                ns("source"), 
                "Data Source", 
                width = 4, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            dinucRegionsSelectorBoxUI(
                ns("dinucRegions")
            ),
            NULL
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("correlationPlotBox"), 
                "Correlation Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns("umapPlotBox"), 
                "UMAP Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            clusterProfilePlotBoxUI(
                ns("clusterProfilePlotBox"), 
                width = 4
            ),
            NULL
        ),
        regionExpansionUI(
            ns("regionExpansion")
        ),
        NULL
    )
}
