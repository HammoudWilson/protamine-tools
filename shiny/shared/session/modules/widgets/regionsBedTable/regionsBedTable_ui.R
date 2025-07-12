#----------------------------------------------------------------------
# UI components for the regionsBedTable widget module
#----------------------------------------------------------------------

# module ui function
regionsBedTableUI <- function(id, width = 4) {
    
    # initialize namespace
    ns <- NS(id)
    
    # box with the table
    column(
        width = width,
        style = "margin: 0; padding: 0;",
        bufferedTableUI(
            ns("table"),
            title = "Enrichment BED Files",
            width = 12,
            status = 'primary',
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = FALSE
        ),
        box(
            width = 12,
            status = "primary",
            solidHeader = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            uiOutput(ns("regionsSummary"))
        )
    )
}
