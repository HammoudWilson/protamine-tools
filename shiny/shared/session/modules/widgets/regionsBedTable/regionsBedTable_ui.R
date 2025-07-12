#----------------------------------------------------------------------
# UI components for the regionsBedTable widget module
#----------------------------------------------------------------------

# module ui function
regionsBedTableUI <- function(id, width = 4, exposePartialOverlap = FALSE) {
    
    # initialize namespace
    ns <- NS(id)
    
    # box with the table
    column(
        width = width,
        style = "margin: 0; padding: 0;",
        if(exposePartialOverlap) box(
            width = 12,
            status = "primary",
            solidHeader = FALSE,
            collapsible = FALSE,
            style = "padding-left: 10px; padding-right: 10px;",
            checkboxInput(
                ns("allow_partial_overlap"), 
                "Allow Partial Overlap", 
                value = TRUE
            )
        ) else NULL,
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
            style = "padding-left: 10px; padding-right: 10px;",
            uiOutput(ns("regionsSummary"))
        )
    )
}
