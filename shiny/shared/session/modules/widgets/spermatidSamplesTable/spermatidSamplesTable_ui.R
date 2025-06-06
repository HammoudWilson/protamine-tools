#----------------------------------------------------------------------
# UI components for the spermatidSamplesTable widget module
#----------------------------------------------------------------------

# module ui function
spermatidSamplesTableUI <- function(id, width = 4) {
    
    # initialize namespace
    ns <- NS(id)

    # known stages
    stages <- c("mESC", "early_RS", "int_RS", "late_RS",
                "earliest_ES", "early_ES", "int_ES", "late_ES")

    tagList(
        # box with table filter inputs
        box(
            width = width,
            title = NULL,
            status = "primary",
            solidHeader = FALSE,
            collapsible = TRUE,
            column(
                width = 2,
                style = "margin-bottom: 5px;",
                textInput(
                    ns("batchFilter"),
                    "Batch",
                    value = ""
                )
            ),
            column(
                width = 10,
                style = "margin-bottom: 5px;",
                checkboxGroupInput(
                    ns("stageFilter"),
                    "Stages",
                    choices = stages,
                    selected = c(),
                    inline = TRUE
                )
            ),
            NULL
        ),

        # box with the table
        box(
            width = width,
            title = "Spermatid ATAC-seq Samples",
            status = 'primary',
            solidHeader = TRUE,
            collapsible = FALSE,
            DTOutput(ns("table"))
        ),

        # statistic summary all samples
        box(
            width = width,
            title = NULL,
            status = 'primary',
            solidHeader = FALSE,
            collapsible = FALSE,
            radioButtons(
                ns("summaryPlotColumn"),
                label = "Summary Plot Column",
                choices = c("primary_count","fraction_spike"),
                selected = "primary_count",
                inline = TRUE
            )
        ),
        staticPlotBoxUI(
            ns("summaryPlot"), 
            "Summary Plot, All Samples",
            width = width,
            status = "primary",
            collapsible = TRUE,
            collapsed = TRUE,
            solidHeader = TRUE
        ),
    )
}
