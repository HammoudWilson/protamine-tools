#----------------------------------------------------------------------
# UI components for the paSampleCounts appStep module
#----------------------------------------------------------------------

# module ui function
paSampleCountsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$paSampleCounts)

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
                width = 12, 
                collapsible = TRUE,
                inFluidRow = FALSE
            ),
            # table of all sample counts by filtering stage
            bufferedTableUI(
                ns("paSampleCountsTable"),
                "Sample Counts by Filtering Stage",
                width = 12,
                status = "primary",
                collapsible = FALSE,
                solidHeader = TRUE,
                downloadable = TRUE
            ),
            NULL
        ),
        NULL
    )
}
