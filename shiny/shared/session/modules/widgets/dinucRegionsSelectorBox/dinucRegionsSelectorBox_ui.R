#----------------------------------------------------------------------
# UI components for the dinucRegionsSelectorBox widget module
#----------------------------------------------------------------------

# module ui function
dinucRegionsSelectorBoxUI <- function(id, width = 4, includeUmap = TRUE) {

    # initialize namespace
    ns <- NS(id)

    # widget UI elements, populate as needed
    # e.g. textInput(ns("xxx"), "XXX")
    # see mdi-apps-framework documentation for useful MDI elements
    box(
        width = width,
        status = "primary",
        collapsible = FALSE,
        title = NULL,
        selectInput(
            ns("index_stage"), 
            "Index Stage", 
            choices = c("overlap_group"),
            selected = "overlap_group",
            width = "100%"
        ),
        div(style = "margin-top: 5px; margin-bottom: 5px;", textInput(
            ns("user_filter"),
            "Where",
            placeholder = "Enter a region filter expression, e.g., early_RS > 5",
            width = "100%"
        )),
        radioButtons(
            ns("rpkm_scaling"), 
            "RPKM Scaling", 
            choices = c("scaled", "unscaled"),
            selected = "scaled",
            inline = TRUE,
            width = "100%"
        ),
        if(includeUmap) div(style = "margin-top: 5px;", radioButtons(
            ns("umap_metric"), 
            "UMAP Metric", 
            choices = c("euclidean", "correlation"),
            selected = "correlation",
            inline = TRUE,
            width = "100%"
        )) else NULL,
        checkboxInput(
            ns("include_unpassed_regions"), 
            "Include Unpassed Regions", 
            value = FALSE
        ),
        checkboxInput(
            ns("allow_partial_overlap"), 
            "Allow Partial Overlap", 
            value = FALSE
        ),
        NULL
    )
}
