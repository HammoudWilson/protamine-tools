#----------------------------------------------------------------------
# UI components for the normalization appStep module
#----------------------------------------------------------------------

# module ui function
normalizationUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalization)

    # UI functions
    plotBox_ <- function(title, ui, collapsed = FALSE) {
        box(
            title = title,
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            collapsed = collapsed,
            ui
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
                spermatidSamplesTableUI(
                    ns("spermatidSamplesTable"),
                    width = 12
                )
            ),
            # (normalized) bin count vs. bin GC content
            column(
                width = 6,
                # box with heatmap setting inputs
                box(
                    width = 12,
                    title = NULL,
                    status = "primary",
                    solidHeader = FALSE,
                    collapsible = TRUE,
                    column(
                        width = 6,
                        style = "margin-bottom: 5px;",
                        checkboxInput(
                            ns("normalizeTn5Site"),
                            "Apply Tn5 Site Normalization",
                            value = FALSE
                        )
                    ),
                    NULL
                ),
                # two matrix plots
                plotBox_(
                    "GC Bias Plot (after any normalization)",
                    interactiveScatterplotUI(ns("gcBiasPlot"), height = '400px'), 
                    collapsed = FALSE
                ),
                plotBox_(
                    "GC Bias Residuals",
                    interactiveScatterplotUI(ns("gcResidualBiasPlot"), height = '400px'),
                    collapsed = FALSE
                ),
                staticPlotBoxUI(
                    ns("tn5vsGCPlot"), 
                    "Tn5 Weighting vs. Bin GC",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                staticPlotBoxUI(
                    ns("wgtDistPlot"), 
                    "Tn5 Weight Distribution",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                NULL
            ),
            NULL
        ),
        NULL
    )
}
