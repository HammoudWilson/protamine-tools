#----------------------------------------------------------------------
# UI components for the dinucRegions appStep module
#----------------------------------------------------------------------

# module ui function
dinucRegionsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$dinucRegions)

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

        # data source selectors
        fluidRow(
            dataSourceTableUI(
                ns("source"), 
                "Data Source", 
                width = 4, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            box(
                width = 4,
                solidHeader = FALSE,
                collapsible = FALSE,
                title = NULL,
                selectInput(
                    ns("index_stage"), 
                    "Index Stage", 
                    choices = c("overlap_group"),
                    selected = "overlap_group",
                    width = "100%"
                ),
                radioButtons(
                    ns("rpkm_scaling"), 
                    "RPKM Scaling", 
                    choices = c("scaled", "unscaled"),
                    selected = "scaled",
                    inline = TRUE,
                    width = "100%"
                ),
                radioButtons(
                    ns("umap_metric"), 
                    "UMAP Metric", 
                    choices = c("euclidean", "correlation"),
                    selected = "correlation",
                    inline = TRUE,
                    width = "100%"
                ),
                checkboxInput(
                    ns("include_unpassed_regions"), 
                    "Include Unpassed Regions", 
                    value = FALSE
                ),
                NULL
            ),
            NULL
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("correlationPlotBox"), 
                "Correlation Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns("umapPlotBox"), 
                "UMAP Plot",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            staticPlotBoxUI(
                ns("clusterProfilePlot"), 
                "Cluster Profiles",
                width = 4,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                solidHeader = TRUE
            ),
            NULL
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("regionProfilePlot"), 
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
        ),
        NULL
    )
}
