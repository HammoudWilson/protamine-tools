#----------------------------------------------------------------------
# UI components for the intervalExpansion widget module
#----------------------------------------------------------------------

# module ui function
intervalExpansionUI <- function(id) {
    ns <- NS(id)
    fluidRow(
        staticPlotBoxUI(
            ns("profilePlotBox"), 
            "Selected Interval Profile",
            width = 4,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE,
            solidHeader = TRUE
        ),
        mdiInteractivePlotBoxUI(
            ns("browserPlotBox"), 
            "Selected Interval Footprint",
            width = 8,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE,
            solidHeader = TRUE
        ),
        NULL
    )
}
