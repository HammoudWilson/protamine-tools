#----------------------------------------------------------------------
# server components for the spermatidStageTable widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
spermatidStageTableServer <-  function(
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
samples <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    paInsertSizes(sourceId)$samples
})
stages <- reactive({
    samples <- samples()
    req(samples)
    samples[order(order)][, .( # order,filename_prefix,sample_name,stage,stage_long,replicate_N,R color Value 
        stage_long = stage_long[1],
        nReplicates = .N
    ), by = .(stage)]
})
selectedSpermatidStages <- reactive({
    rows <- selectedRows()
    req(rows)
    stages()[rows, stage]
})

#----------------------------------------------------------------------
# render the selection table
#----------------------------------------------------------------------
output$table <- renderDT(
    {
        x <- stages()
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
    selectedStages = selectedSpermatidStages,
    selectedSamples = reactive({
        samples()[stage %in% selectedSpermatidStages()]
    })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
