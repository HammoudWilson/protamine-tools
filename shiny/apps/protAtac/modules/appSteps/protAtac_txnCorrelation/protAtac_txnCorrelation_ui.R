#----------------------------------------------------------------------
# UI components for the protAtac_txnCorrelation appStep module
#----------------------------------------------------------------------

# module ui function
protAtac_txnCorrelationUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$protAtac_txnCorrelation)

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
            bufferedTableUI(
                ns("sample"),
                "Sample",
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = FALSE
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
                        "GC Residual Z Score", 
                        "Intermediate Insert Size Fraction",
                        "Protamine Enrichment NRLL"
                    ),
                    choiceValues = c(
                        "gcrz", 
                        "iisf",
                        "nrll"
                    ),
                    selected = "gcrz",
                    inline = TRUE,
                    width = "100%"
                )
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("txnCorrelationPlot"), 
                "Transcription Correlation",
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
