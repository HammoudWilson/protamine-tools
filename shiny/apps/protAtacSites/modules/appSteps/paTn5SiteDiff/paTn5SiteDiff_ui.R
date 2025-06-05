#----------------------------------------------------------------------
# UI components for the paTn5SiteDiff appStep module
#----------------------------------------------------------------------

# module ui function
paTn5SiteDiffUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$paTn5SiteDiff)

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
            dataSourceTableUI(
                ns("dataSourceTable"),
                "Data Source", 
                width = 6, 
                collapsible = TRUE,
                inFluidRow = FALSE
            ),
            staticPlotBoxUI(
                ns("dinucEnrichmentPlot"), 
                "Tn5 Site Enrichment by Dinucleotide",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            NULL
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("tn5DiffPlot_CG"), 
                "Tn5 RS vs. ES Difference, CpG",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            staticPlotBoxUI(
                ns("tn5DiffPlot_GC"), 
                "Tn5 RS vs. ES Difference, GpC",
                width = 6,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            NULL
        ),
        NULL
    )
}
