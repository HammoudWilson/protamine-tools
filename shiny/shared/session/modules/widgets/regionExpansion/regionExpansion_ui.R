#----------------------------------------------------------------------
# UI components for the regionExpansion widget module
#----------------------------------------------------------------------

# module ui function
regionExpansionUI <- function(id) {

    # initialize namespace
    ns <- NS(id)

    # widget UI elements, populate as needed
    # e.g. textInput(ns("xxx"), "XXX")
    # see mdi-apps-framework documentation for useful MDI elements
    fluidRow(
        staticPlotBoxUI(
            ns("regionProfilePlotBox"), 
            "Selected Region Profile",
            width = 4,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE,
            solidHeader = TRUE
        ),
        mdiInteractivePlotBoxUI(
            ns("regionPlotBox"), 
            "Selected Region Footprint",
            width = 8,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE,
            solidHeader = TRUE
        ),
        NULL
    )
}
