#----------------------------------------------------------------------
# UI components for the paTn5Sites appStep module
#----------------------------------------------------------------------

# module ui function
paTn5SitesUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$paTn5Sites)

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
        # settings = TRUE,

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
                # box with settings inputs
                box(
                    width = 12,
                    title = NULL,
                    status = "primary",
                    solidHeader = FALSE,
                    collapsible = TRUE,
                    column(
                        width = 12,
                        style = "margin-bottom: 5px;",
                        radioButtons(
                            ns("tn5Site_rankBy"),
                            "Rank Tn5 Sites By",
                            choices = c("frequency", "enrichment"),
                            selected = "frequency",
                            inline = TRUE
                        )
                    )
                ),
                staticPlotBoxUI(
                    ns("tn5CdfPlot"), 
                    "Tn5 CDF Plot",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = FALSE,
                    solidHeader = TRUE
                ),
                bufferedTableUI(
                    ns("tn5SiteTable"),
                    "Ranked Tn5 Sites",
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
