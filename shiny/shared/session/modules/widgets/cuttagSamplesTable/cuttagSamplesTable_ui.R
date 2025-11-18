#----------------------------------------------------------------------
# UI components for the cuttagSamplesTable widget module
#----------------------------------------------------------------------

# module ui function
cuttagSamplesTableUI <- function(id, width = 4) {
    
    # initialize namespace
    ns <- NS(id)

    # box with the table
    box(
        width = width,
        title = "Spermatid Cut&Tag Samples",
        status = 'primary',
        solidHeader = TRUE,
        collapsible = FALSE,
        DTOutput(ns("table"))
    )
}
