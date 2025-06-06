#----------------------------------------------------------------------
# server components for the spermatidSamplesTable widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
spermatidSamplesTableServer <-  function(
    id, 
    samples, # a reactive that returns a samples metadata data.table
    selection = "multiple",
    n_ref_wgt_is_gc_smp = NULL # a reactive that returns similarly named object from collate output, when available
) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# track the table, i.e., data source selection
#----------------------------------------------------------------------
selectedRows <- rowSelectionObserver('table', input)
allSamples <- reactive({
    samples <- samples()
    req(samples)
    samples[order(staging_order)]
})
filteredSamples <- reactive({
    filteredSamples <- allSamples()
    req(filteredSamples)
    if (input$batchFilter != "") {
        batch_ <- tryCatch({
            as.integer(trimws(input$batchFilter))
        }, error = function(e) {
            NA_integer_
        })
        if (!is.na(batch_)) {
            filteredSamples <- filteredSamples[batch == batch_]
        }
    }
    if (length(input$stageFilter) > 0) {
        filteredSamples <- filteredSamples[stage %in% input$stageFilter]
    }
    filteredSamples
})
selectedSpermatidSamples <- reactive({
    rows <- selectedRows()
    req(rows)
    filteredSamples()[rows][order(staging_order)]
})

#----------------------------------------------------------------------
# render the selection table
#----------------------------------------------------------------------
getNInserts <- function(refType, sampleNames){
    if (is.null(n_ref_wgt_is_gc_smp)) {
        NA_integer_
    } else {
        n_is_gc_smp <- n_ref_wgt_is_gc_smp()[[refType]]$observed
        sapply(sampleNames, function(sample_name) {
            sum(n_is_gc_smp[,,sample_name], na.rm = TRUE)
        })
    }
}
tableData <- reactive({
    x <- filteredSamples()[,
        .(
            batch   = batch, 
            prefix  = filename_prefix, 
            name    = sample_name, 
            staging = staging, 
            stage   = stage,
            primary  = getNInserts("primary",  sample_name),
            spike_in = getNInserts("spike_in", sample_name)
        )
    ]
    x[, frac_spike := round(spike_in / (primary + spike_in), 3)]
    stopSpinner(session)
    x
})
output$table <- renderDT(
    { tableData() },
    options = list(
        paging = FALSE,
        searching = FALSE  
    ),
    class = "display table-compact-4",
    selection = selection,
    editable = FALSE, 
    rownames = FALSE # must be true for editing to work, not sure why (datatables peculiarity)
)

#----------------------------------------------------------------------
# render the all summary plot
#----------------------------------------------------------------------
summaryPlot <- staticPlotBoxServer(
    "summaryPlot",
    maxHeight = "400px",
    create = function() {
        samples <- tableData()
        req(samples)
        columns <- list(
            primary_count = "primary",
            fraction_spike = "frac_spike"
        )
        d <- samples[[columns[[input$summaryPlotColumn]]]]
        if(input$summaryPlotColumn == "primary_count") d <- d / 1e6 # convert to millions
        ylab <- if(input$summaryPlotColumn == "primary_count") "Primary Inserts (M)" 
                else "Fraction Spike-in"
        names(d) <- samples$staging
        par(mar = c(5, 4, 0, 0) + 0.1)
        barplot(d, las=2, ylab = ylab, cex.names = 0.75)
    }

)

#----------------------------------------------------------------------
# return a reactive populated with the id(s) of the selected source(s)
#----------------------------------------------------------------------
list(
    allSamples = allSamples,
    selectedSamples = selectedSpermatidSamples
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
