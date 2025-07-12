#----------------------------------------------------------------------
# UI components for the clusterProfilePlotBox widget module
#----------------------------------------------------------------------

# module ui function
clusterProfilePlotBoxUI <- function(id, width = 4) {
    ns <- NS(id)
    staticPlotBoxUI(
        ns("plot"), 
        "Cluster Profiles",
        width = width,
        status = "primary",
        collapsible = TRUE,
        collapsed = FALSE,
        solidHeader = TRUE
    )
}
