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
            width = 6,
            radioButtons(
                ns(id),
                paste0("Score Type, ", axis, "-axis"),
                choiceNames  = c(
                    "Bin GC", 
                    "Pro-seq Txn",
                    "HiC Compartment",
                    "ATAC Stage Mean", 
                    "ATAC GC Resid Z",
                    "Insert Size NRLL",
                    "H2B",
                    "H4",
                    "H4ac",
                    "H3K27me3"
                ),
                choiceValues = c(
                    "gc", 
                    "txn", 
                    "hic",
                    "stgm",
                    "gcrz_obs", 
                    "nrll",
                    "H2B",
                    "H4",
                    "H4ac",
                    "H3K27me3"
                ),
                selected = selected,
                inline = TRUE,
                width = "100%"
            )
        )
    }
    scoreStageSelecton <- function(id, axis){
        box(
            title = NULL,
            status = "primary",
            solidHeader = FALSE,
            collapsible = FALSE,
            width = 2,
            selectInput(
                ns(id),
                paste0("Spermatid Stage, ", axis, "-axis"),
                choices = c(
                    "early_RS",
                    "int_RS",
                    "late_RS",
                    "earliest_ES",
                    "early_ES",
                    "int_ES",
                    "late_ES",
                    "round - elong"
                ),
                selected = "early_RS",
                multiple = FALSE,
                width = "100%"
            )
        )
    }
    scoreCheckboxGroup <- function(id){
        box(
            title = NULL,
            status = "primary",
            solidHeader = FALSE,
            collapsible = FALSE,
            width = 6,
            checkboxGroupInput(
                ns(id),
                "PCA Score Types",
                choiceNames  = c(
                    "Bin GC", 
                    "Pro-seq Txn",
                    "HiC Compartment",
                    "ATAC Stage Mean", 
                    "ATAC GC Resid Z",
                    "H2B",
                    "H4",
                    "H4ac",
                    "H3K27me3"
                ),
                choiceValues = c(
                    "gc_z", 
                    "txn", 
                    "hic",
                    "stgm",
                    "gcrz_obs", 
                    "H2B",
                    "H4",
                    "H4ac",
                    "H3K27me3"
                ),
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

        # appStep UI elements, populate as 
        fluidRow(
            dataSourceTableUI(
                ns("dataSourceTable"),
                "Data Source", 
                width = 6, 
                collapsible = TRUE,
                inFluidRow = FALSE
            )
        ),
        fluidRow(
            scoreRadioButtons("xScoreType", "X", "hic"),
            scoreStageSelecton("xScoreStage", "X")
        ),
        fluidRow(
            scoreRadioButtons("yScoreType", "Y", "gc"),
            scoreStageSelecton("yScoreStage", "Y")
        ),
        fluidRow(
            scoreRadioButtons("zScoreType", "Z", "txn"),
            scoreStageSelecton("zScoreStage", "Z")
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("correlationPlot"), 
                "Score Correlation",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = FALSE
            )
        ),
        fluidRow(
            scoreCheckboxGroup("pcaScoreTypes")
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("pcaPlot"), 
                "PCA Plot",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                solidHeader = TRUE,
                collapsed = TRUE
            )
        ),
        NULL
    )
}
