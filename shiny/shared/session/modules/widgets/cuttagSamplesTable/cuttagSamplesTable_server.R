#----------------------------------------------------------------------
# server components for the cuttagSamplesTable widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
cuttagSamplesTableServer <-  function(
    id, 
    samples, # a reactive that returns a samples metadata data.table
    selection = "multiple"
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
    samples[order(class_order, target_order, staging_order)]
})
selectedSamples <- reactive({
    rows <- selectedRows()
    req(rows)
    allSamples()[rows][order(class_order, target_order, staging_order)]
})

#----------------------------------------------------------------------
# render the selection table
#----------------------------------------------------------------------
tableData <- reactive({
    x <- allSamples()[,
        .( 
            name    = sample_name, 
            antibody_target = antibody_target,
            staging = staging, 
            stage   = stage
        )
    ]
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
# return a reactive populated with the id(s) of the selected source(s)
#----------------------------------------------------------------------
list(
    allSamples = allSamples,
    selectedSamples = selectedSamples
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
