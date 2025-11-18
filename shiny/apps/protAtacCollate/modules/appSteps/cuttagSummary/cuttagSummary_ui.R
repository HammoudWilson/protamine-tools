#----------------------------------------------------------------------
# UI components for the cuttagSummary appStep module
#----------------------------------------------------------------------

# module ui function
cuttagSummaryUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$cuttagSummary)

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
        settings = FALSE,

        # appStep UI elements, populate as needed
        fluidRow(
            # data selection inputs
            column(
                width = 6,
                dataSourceTableUI(
                    ns("dataSourceTable"),
                    "Data Source", 
                    width = 12, 
                    collapsible = TRUE,
                    inFluidRow = FALSE
                ),
                cuttagSamplesTableUI(
                    ns("cuttagSamplesTable"),
                    width = 12
                )
            ),
            column(
                width = 6,
                box(
                    title = NULL,
                    status = "primary",
                    solidHeader = FALSE,
                    collapsible = FALSE,
                    width = 12,
                    radioButtons(
                        ns("scoreType"),
                        "Score Type",
                        choiceNames  = c(
                            "H2B",
                            "H4",
                            "H3K27me3",
                            "H4ac"
                        ),
                        choiceValues = c(
                            "H2B",
                            "H4",
                            "H3K27me3",
                            "H4ac"
                        ),
                        selected = "H2B",
                        inline = TRUE,
                        width = "100%"
                    )
                ),
                staticPlotBoxUI(
                    ns("sampleDistributionPlot"), 
                    "Sample Distributions",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = FALSE
                ),
                staticPlotBoxUI(
                    ns("stageDistributionPlot"), 
                    "Stage Distributions",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = FALSE
                )
            ),
            NULL
        ),
        NULL
    )
}
