#----------------------------------------------------------------------
# UI components for the spermatidStageTable widget module
#----------------------------------------------------------------------

# module ui function
spermatidStageTableUI <- function(id, width = 4) {
    
    # initialize namespace
    ns <- NS(id)
    
    # box with the table
    box(
        width = width,
        title = "Spermatid Stages",
        status = 'primary',
        solidHeader = TRUE,
        collapsible = FALSE,
        DTOutput(ns("table"))
    )
}
