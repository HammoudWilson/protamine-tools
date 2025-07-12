#----------------------------------------------------------------------
# UI components for the dinucChainsSelectorBox widget module
#----------------------------------------------------------------------

# module ui function
dinucChainsSelectorBoxUI <- function(id, width = 4, includeUmap = TRUE) {
    ns <- NS(id)
    box(
        width = width,
        status = "primary",
        collapsible = FALSE,
        title = NULL,
        style = "padding-left: 10px; padding-right: 10px;",
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
            placeholder = "Enter a filter expression, e.g., early_RS > 5",
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
            ns("include_unpassed"), 
            "Include Unpassed Intervals", 
            value = FALSE
        ),
        NULL
    )
}
