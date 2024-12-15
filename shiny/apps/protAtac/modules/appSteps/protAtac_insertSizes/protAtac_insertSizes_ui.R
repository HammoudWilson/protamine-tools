#----------------------------------------------------------------------
# UI components for the protAtac_insertSizes appStep module
#----------------------------------------------------------------------

# module ui function
protAtac_insertSizesUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$protAtac_insertSizes)

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
            dataSourceTableUI(
                ns("dataSourceTable"),
                "Data Source", 
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            spermatidStageTableUI(
                ns("spermatidStageTable"),
                width = 6
            )
        ),
        fluidRow(
            box(
                title = "Plot Options",
                status = "primary",
                solidHeader = TRUE,
                collapsible = FALSE,
                width = 12,
                column(
                    width = 2,
                    checkboxInput(
                        ns("aggregateSamplesByStage"),
                        "Aggregate Samples By Stage",
                        value = TRUE
                    )
                ),
                column(
                    width = 10,
                    checkboxInput(
                        ns("normalizeToSpikeIn"),
                        "Normalize To Spike In Samples By Stage",
                        value = FALSE
                    )
                )
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("insertSizesPlot_genome"), 
                "Spermatid Insert Sizes",
                width = 6,
                status = "primary",
                collapsible = FALSE,
                solidHeader = TRUE
            ),
            staticPlotBoxUI(
                ns("insertSizesPlot_spike_in"), 
                "Spike-in Insert Sizes",
                width = 6,
                status = "primary",
                collapsible = FALSE,
                solidHeader = TRUE
            )
        ),
        NULL
    )
}
