#----------------------------------------------------------------------
# UI components for the scoreCorrelation appStep module
#----------------------------------------------------------------------

# module ui function
scoreCorrelationUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$scoreCorrelation)

    # score selectors
    scoreRadioButtons <- function(id, axis, selected){
        box(
            title = NULL,
            status = "primary",
            solidHeader = FALSE,
            collapsible = FALSE,
            width = 12,
            radioButtons(
                ns(id),
                paste0("Score Type, ", axis, "-axis"),
                choiceNames  = c(
                    "Bin GC", 
                    "Pro-seq Txn",
                    "Bin Stage Mean", 
                    "HiC Compartment Score",
                    "GC Resid Z Obs Delta",
                    "Protamine NRLL Delta"
                ),
                choiceValues = c(
                    "gc", 
                    "txn",
                    "stgm", 
                    "hic",
                    "gcrz_obs", 
                    "nrll"
                ),
                selected = selected,
                inline = TRUE,
                width = "100%"
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
        settings = TRUE,

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
                scoreRadioButtons("xScoreType", "X", "hic"),
                scoreRadioButtons("yScoreType", "Y", "gc"),
                scoreRadioButtons("zScoreType", "Z", "txn"),
                staticPlotBoxUI(
                    ns("correlationPlot"), 
                    "Score Correlation",
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
