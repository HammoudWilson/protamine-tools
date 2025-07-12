#----------------------------------------------------------------------
# UI components for the regionEnrichment appStep module
#----------------------------------------------------------------------

# module ui function
regionEnrichmentUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$regionEnrichment)

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
            # BED file upload
            column(
                width = 6,
                style = "margin: 0; padding: 0;",
                regionsBedTableUI(
                    ns("regionsBedTable"),
                    width = 12
                ),
                box(
                    title = NULL,
                    status = "primary",
                    solidHeader = FALSE,
                    collapsible = FALSE,
                    width = 12,
                    radioButtons(
                        ns("enrichmentScoreType"),
                        "Enrichment Score Type",
                        choices = c(
                            "gcrz_obs",
                            "gcrz_wgt",
                            "nrll"
                        ),
                        selected = "gcrz_obs",
                        inline = TRUE
                    )
                ),
                staticPlotBoxUI(
                    ns("bySamplePlot"), 
                    "Sample Enrichment",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = TRUE
                ),
                staticPlotBoxUI(
                    ns("byStagePlot"), 
                    "Stage Enrichment",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = TRUE
                ),
                staticPlotBoxUI(
                    ns("byStageTypePlot"), 
                    "Stage Type Enrichment",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE,
                    collapsed = TRUE
                ),
                NULL
            ),
            NULL
        ),
        NULL
    )
}
