#----------------------------------------------------------------------
# server components for the spermatidSamplesTable widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
spermatidSamplesTableServer <-  function(
    id, 
    sourceId, # a reactive that returns the id of one selected source
    selection = "multiple"
) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# track the table, i.e., data source selection
#----------------------------------------------------------------------
selectedRows <- rowSelectionObserver('table', input)
allSamples <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paInsertSizes(sourceId)$samples[order(staging_order)]
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
    insertSizes <- paInsertSizes(sourceId())$insertSizes
    sapply(sampleNames, function(x) {
        sum(insertSizes[[x]][[refType]], na.rm = TRUE)
    })
}

output$table <- renderDT(
    {
        x <- filteredSamples()[,
            .(
                batch   = batch, 
                prefix  = filename_prefix, 
                name    = sample_name, 
                staging = staging, 
                stage   = stage,
                genome   = getNInserts("genome",   sample_name),
                spike_in = getNInserts("spike_in", sample_name)
            )
        ]
        x[, frac_spike := round(spike_in / (genome + spike_in), 3)]
        stopSpinner(session)
        x
    },
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
