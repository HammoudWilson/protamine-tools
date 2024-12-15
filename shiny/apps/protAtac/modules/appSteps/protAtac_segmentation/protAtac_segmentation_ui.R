#----------------------------------------------------------------------
# UI components for the protAtac_segmentation appStep module
#----------------------------------------------------------------------

# module ui function
protAtac_segmentationUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$protAtac_segmentation)

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
            dataSourceTableUI(
                ns("source"),
                "Data Source", 
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            bufferedTableUI(
                ns("trainingRegions"),
                "Training Regions",
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = FALSE, 
                downloadable = TRUE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("gcrzDistributions"),
                "GCRZ Distributions",
                width = 4,
                collapsible = FALSE
            ),
            staticPlotBoxUI(
                ns("iisfDistributions"),
                "IISF Distributions",
                width = 4,
                collapsible = FALSE
            ),
            staticPlotBoxUI(
                ns("nrllDistributions"),
                "NRLL Distributions",
                width = 4,
                collapsible = FALSE
            )
        ),
        fluidRow(
            column(
                width = 3, 
                bsButton(
                    ns("clearTrainingRegions"),
                    "Clear Training Regions",
                    style = "danger"
                )
            ),
            column(
                width = 3,
                bsButton(
                    ns("solveHMM"),
                    "Run HMM Segmentation",
                    style = "success"
                )
            ),
            column(
                width = 3,
                textInput(
                    ns("trainingSetName"),
                    "Training Set Name",
                    value = ""
                )
            )
        ),
        NULL
    )
}
