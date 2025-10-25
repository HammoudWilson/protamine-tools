#----------------------------------------------------------------------
# UI components for the lagAnalysis appStep module
#----------------------------------------------------------------------

# module ui function
lagAnalysisUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$lagAnalysis)

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
            NULL
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("variogramPlotBox"), 
                "Variogram",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns("autocorrelationPlotBox"), 
                "Autocorrelation (sort of)",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns("diffPlotBox"), 
                "Stage Mean Difference",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            NULL
        ),
        NULL
    )
}
