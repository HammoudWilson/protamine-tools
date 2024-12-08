#----------------------------------------------------------------------
# UI components for the protAtac_normalizeGC appStep module
#----------------------------------------------------------------------

# module ui function
protAtac_normalizeGCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$protAtac_normalizeGC)

    # UI functions
    plotBox_ <- function(title, column1UI, column2UI = NULL, column3UI = NULL){
        box(
            title = title,
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            # collapsible = TRUE,
            column(
                width = 3,
                column1UI
            ),
            column(
                width = 9,
                column2UI
            ),
            column(
                width = 12,
                column3UI
            )
        )
    }

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
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            bufferedTableUI(
                ns("sample"),
                "Sample",
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = FALSE
            )
        ),

        # GC bias plots
        fluidRow(
            plotBox_(
                "GC Bias",
                tags$div(
                    numericInput(ns("nBiasBins"), "# Plotted Bins", value = 10000, min = 1000, max = 50000, step = 1000)
                ),
                interactiveScatterplotUI(ns("gcBiasPlot"), height = '400px')
            ),
            plotBox_(
                title = "GC Bias Residuals",
                tags$div(
                    numericInput(ns("nResidualBiasBins"), "# Plotted Bins", value = 10000, min = 1000, max = 50000, step = 1000)
                ),
                interactiveScatterplotUI(ns("gcResidualBiasPlot"), height = '400px')
            )
        ),
        NULL
    )
}
