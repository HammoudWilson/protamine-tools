#----------------------------------------------------------------------
# UI components for the modelParameters appStep module
#----------------------------------------------------------------------

# module ui function
modelParametersUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$modelParameters)

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

        # appStep UI elements, populate as needed
        fluidRow(
            # data selection inputs
            column(
                width = 6,
                dataSourceTableUI(
                    ns("dataSourceTable"),
                    "Data Source", 
                    width = 12, 
                    collapsible = TRUE,
                    inFluidRow = FALSE
                ),
                spermatidSamplesTableUI(
                    ns("spermatidSamplesTable"),
                    width = 12
                )
            ),
            # insert size by GC percent heatmaps
            column(
                width = 6,
                # # box with heatmap setting inputs
                # box(
                #     width = 12,
                #     title = NULL,
                #     status = "primary",
                #     solidHeader = FALSE,
                #     collapsible = TRUE,
                #     column(
                #         width = 3,
                #         style = "margin-bottom: 5px;",
                #         checkboxInput(
                #             ns("scaleHeatmapByInsertSize"),
                #             "Scale By Size",
                #             value = FALSE
                #         )
                #     ),
                #     column(
                #         width = 3,
                #         style = "margin-bottom: 5px;",
                #         checkboxInput(
                #             ns("nucleosomesOnly"),
                #             "Nucleosomes Only",
                #             value = FALSE
                #         )
                #     ),
                #     NULL
                # ),
                # mappability plot
                staticPlotBoxUI(
                    ns("mappabilityPlot"), 
                    "Mappability",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                # mappability distribution plot
                staticPlotBoxUI(
                    ns("mappabilityDistributionPlot"), 
                    "Mappability Distribution",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                # GC plot
                staticPlotBoxUI(
                    ns("gcPlot"), 
                    "GC Plot",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                # GC plot
                staticPlotBoxUI(
                    ns("gcPlot2"), 
                    "GC Plot / Mappability",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                # Tn5 plots
                staticPlotBoxUI(
                    ns("tn5Plot"), 
                    "Tn5 Site Plot",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                staticPlotBoxUI(
                    ns("tn5CdfPlot"), 
                    "Tn5 CDF Plot",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE
                ),
                # staticPlotBoxUI(
                #     ns("tn5DiffPlot"), 
                #     "Tn5 RS vs. ES Diff",
                #     width = 12,
                #     status = "primary",
                #     collapsible = TRUE,
                #     solidHeader = TRUE
                # ),
                box(
                    width = 12,
                    title = "Tn5 RS vs. ES Diff", # arguments passed to shinydashboard::box()
                    status = "primary",
                    solidHeader = TRUE, 
                    collapsible = TRUE, 
                    collapsed = FALSE,
                    interactiveScatterplotUI(ns("tn5DiffPlot"))
                ),
                bufferedTableUI(
                    ns("tn5Table"),
                    "Tn5 Sites",
                    width = 12,
                    status = "primary",
                    collapsible = TRUE,
                    solidHeader = TRUE
                ),
                # staticPlotBoxUI(
                #     ns("insertSizeMatrix_spike_in"), 
                #     "Spike-in Insert Size vs. Insert GC",
                #     width = 12,
                #     status = "primary",
                #     collapsible = TRUE,
                #     collapsed = TRUE,
                #     solidHeader = TRUE
                # ),
                # box(
                #     width = 12,
                #     title = NULL,
                #     status = "primary",
                #     solidHeader = FALSE,
                #     collapsible = TRUE,
                #     column(
                #         width = 3,
                #         checkboxInput(
                #             ns("aggregateByStage"),
                #             "Aggregate By Stage",
                #             value = FALSE
                #         )
                #     ),
                #     column(
                #         width = 4,
                #         checkboxInput(
                #             ns("normalizeToSpikeIn"),
                #             "Normalize To Spike In",
                #             value = FALSE
                #         )
                #     ),
                #     NULL
                # ),
                # # two insert size distribution plots
                # staticPlotBoxUI(
                #     ns("insertSizesPlot_primary"), 
                #     "Spermatid Insert Sizes",
                #     width = 12,
                #     status = "primary",
                #     collapsible = TRUE,
                #     solidHeader = TRUE
                # ),
                # staticPlotBoxUI(
                #     ns("insertSizesPlot_spike_in"), 
                #     "Spike_In Insert Sizes",
                #     width = 12,
                #     status = "primary",
                #     collapsible = TRUE,
                #     solidHeader = TRUE
                # ),
                NULL
            ),
            NULL
        ),
        NULL
    )
}
