#----------------------------------------------------------------------
# UI components for the vPlots appStep module
#----------------------------------------------------------------------

# module ui function
vPlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$vPlots)

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
        settings = TRUE,

        # data source selectors
        fluidRow(
            dataSourceTableUI(
                ns("source"), 
                "Data Source", 
                width = 4, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            bufferedTableUI(
                ns("samples1"),
                "Sample(s) #1",
                width = 4,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE
            ),
            bufferedTableUI(
                ns("samples2"),
                "Sample(s) #2",
                width = 4,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE
            )
        ),

        # V plots
        fluidRow(
            staticPlotBoxUI(
                ns("vPlot1"), 
                "Sample(s) #1 V Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("vPlot2"), 
                "Sample(s) #2 V Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("vPlot_diff"), 
                "Sample(s) #2 - Sample(s) #1 V Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            )
        ),

        NULL
    )
}
