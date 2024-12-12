#----------------------------------------------------------------------
# UI components for the protAtac_scoreSummary appStep module
#----------------------------------------------------------------------

# module ui function
protAtac_scoreSummaryUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$protAtac_scoreSummary)

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
                title = "Score Type Switchboard",
                status = "primary",
                solidHeader = TRUE,
                collapsible = FALSE,
                width = 12,
                radioButtons(
                    ns("scoreType"),
                    "Score Type",
                    choiceNames  = c(
                        "Fraction GC", 
                        "Transcription CPM", 
                        "Counts Per Million", 
                        "GC-Residual Z-Score", 
                        "Subnucleosomal Insert Fraction", 
                        "Subnucleosomal vs. Nucleosomal NRLL"
                    ),
                    choiceValues = c(
                        "gc", 
                        "txn", 
                        "cpm", 
                        "gcrz", 
                        "snif", 
                        "nrll"
                    ),
                    selected = "gc",
                    inline = TRUE,
                    width = "100%"
                )
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("sampleDistributionPlot"), 
                "Sample Distributions",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("stageDistributionPlot"), 
                "Stage Distributions",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("stageTypeDistributionPlot"), 
                "Stage Type Distributions",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("deltaDistributionPlot"), 
                "Stage Type Delta Distribution",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            )
        ),
        NULL
    )
}
