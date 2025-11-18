#----------------------------------------------------------------------
# UI components for the downloader appStep module
#----------------------------------------------------------------------

# module ui function
downloaderUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$downloader)

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
            style = "font-weight: bold; margin: 20px;",
            column(
                width = 12,
                tags$p("Select a Data Source."),
                tags$p("Use the options to configure the kind of table to create."),
                tags$p("Provide a path to an existing directory on the server where data will be written."),
                tags$p("Click the button to write the table."),
                tags$p("BE PATIENT: these are BIG tables and it will take several minutes to create and write them.")
            )
        ),
        fluidRow(
            column(
                width = 5,
                dataSourceTableUI(
                    ns("source"), 
                    "Data Source", 
                    width = 12, 
                    collapsible = TRUE,
                    inFluidRow = FALSE
                )
            ),
            column(
                width = 7,
                checkboxInput(
                    ns("includeHeader"),
                    label = "Include Header",
                    value = TRUE
                ),
                selectInput(
                    ns("fileFormat"),
                    label = "File Format",
                    choices = c("bed.gz", "rds"),
                    selected = "bed.gz"
                ),
                selectInput(
                    ns("scoreType"),
                    label = "Score Type",
                    choices = c("gcrz_obs", "gcrz_wgt", "nrll", "dinuc_chains","H2B","H4","H4ac","H3K27me3"),
                    selected = "gcrz_obs"
                ),
                textInput(
                    ns("outputDir"),
                    label = "Output Directory",
                    value = "", 
                    width = "100%"
                ),
                tags$div(
                    style = "margin-top: 20px;",
                    actionButton(
                        ns("writeData"), 
                        label = "Write Data to File")
                )
            ),
            NULL
        ),
        NULL
    )
}
