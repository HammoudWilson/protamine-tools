#----------------------------------------------------------------------
# UI components for the insertSizes appStep module
#----------------------------------------------------------------------

# module ui function
insertSizesUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$insertSizes)

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
                spermatidSamplesTableUI(
                    ns("spermatidSamplesTable"),
                    width = 12
                )
            ),
            # insert size by GC percent heatmaps
            column(
                width = 6,
                box(
                    width = 12,
                    title = NULL,
                    status = "primary",
                    solidHeader = FALSE,
                    collapsible = TRUE,
                    column(
                        width = 3,
                        checkboxInput(
                            ns("aggregateByStage"),
                            "Aggregate By Stage",
                            value = FALSE
                        )
                    ),
                    NULL
                ),
                # two insert size distribution plots
                staticPlotBoxUI(
                    ns("insertSizesPlot_peak"), 
                    "Spermatid Peak Insert Sizes",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE
                ),
                staticPlotBoxUI(
                    ns("insertSizesPlot_non_peak"), 
                    "Spermatid Non-Peak Insert Sizes",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE
                ),
                NULL
            ),
            NULL
        ),
        NULL
    )
}
