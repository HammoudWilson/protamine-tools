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
                width = 8, 
                collapsible = TRUE,
                inFluidRow = FALSE
            ),
            NULL
        ),
        fluidRow(
            column(
                width = 4,
                style = "margin: 0; padding: 0;",
                mdiInteractivePlotBoxUI(
                    ns("correlationPlotBox"), 
                    "Correlation Plot",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = TRUE
                ),
                mdiInteractivePlotBoxUI(
                    ns("umapPlotBox"), 
                    "UMAP Plot",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = TRUE
                ),
                staticPlotBoxUI(
                    ns("regionProfilePlot"), 
                    "Selected Region Profile",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = TRUE
                ),
                staticPlotBoxUI(
                    ns("clusterProfilePlot"), 
                    "Cluster Profiles",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = FALSE
                ),
                NULL
            ),
            column(
                width = 8,
                style = "margin: 0; padding: 0;",
                mdiInteractivePlotBoxUI(
                    ns("regionPlotBox"), 
                    "Selected Region Footprint",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = FALSE
                ),
                NULL
            )
        ),
        NULL
    )
}
