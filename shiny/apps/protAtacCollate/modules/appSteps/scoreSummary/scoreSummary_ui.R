#----------------------------------------------------------------------
# UI components for the scoreSummary appStep module
#----------------------------------------------------------------------

# module ui function
scoreSummaryUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$scoreSummary)

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
                    title = NULL,
                    status = "primary",
                    solidHeader = FALSE,
                    collapsible = FALSE,
                    width = 12,
                    radioButtons(
                        ns("scoreType"),
                        "Score Type",
                        choiceNames  = c(
                            "Bin GC", 
                            "Pro-seq Txn", 
                            "Stage Mean", 
                            "HiC Compartment",
                            "GC Resid Z Obs",
                            "GC Resid Z Wgt", 
                            "Protamine NRLL"
                        ),
                        choiceValues = c(
                            "gc", 
                            "txn", 
                            "stgm", 
                            "hic",
                            "gcrz_obs", 
                            "gcrz_wgt", 
                            "nrll"
                        ),
                        selected = "gc",
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
                ),
                staticPlotBoxUI(
                    ns("stageTypeDistributionPlot"), 
                    "Stage Type Distributions",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = FALSE
                ),
                staticPlotBoxUI(
                    ns("deltaDistributionPlot"), 
                    "Stage Type Delta Distribution",
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
