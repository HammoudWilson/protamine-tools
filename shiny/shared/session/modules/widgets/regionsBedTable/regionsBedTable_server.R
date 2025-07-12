#----------------------------------------------------------------------
# server components for the regionsBedTable widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
regionsBedTableServer <-  function(
    id, 
    sourceId # a reactive that returns the id of one selected source
) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# available BED region files
#----------------------------------------------------------------------
regionsBedDir <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    file.path(
        getSourcePackageOption(sourceId, "output", "output-dir"),
        getSourcePackageOption(sourceId, "output", "data-name"),
        "regions_bed"
    )
})
bedFiles <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    regionsBedDir <- regionsBedDir()
    req(regionsBedDir)
    regionClasses <- list.dirs(
        regionsBedDir, 
        recursive = FALSE, 
        full.names = FALSE
    )
    req(length(regionClasses) > 0)
    do.call(rbind, lapply(regionClasses, function(regionClass) {
        regionTypes <- list.files(
            file.path(regionsBedDir, regionClass),
            pattern = "\\.bed$",
            recursive = FALSE, 
            full.names = FALSE
        )
        if (length(regionTypes) == 0) return(NULL)
        data.table(
            regionClass = regionClass,
            regionType  = sub(".bed", "", regionTypes),
            path        = file.path(regionsBedDir, regionClass, regionTypes),
            name        = paste(regionClass, regionTypes, sep = "/")
        )[order(regionClass, regionType)]
    }))
})

#----------------------------------------------------------------------
# render the selection table
#----------------------------------------------------------------------
tableData <- reactive({
    bedFiles <- bedFiles()
    req(bedFiles)
    bedFiles[, .(regionClass, regionType)]
})
table <- bufferedTableServer(
    "table",
    id,
    input,
    tableData,
    selection = "single",
    options = list()
)

#----------------------------------------------------------------------
# collect the selected regions BED
#----------------------------------------------------------------------
bed3Cols <- c("chrom", "start0", "end1")
bed4Cols <- c(bed3Cols, "name")
selectedBed <- reactive({
    selectedRow <- table$rows_selected()
    req(selectedRow)
    bedFiles()[selectedRow]
})
selectedBedData <- reactive({
    selectedBed <- selectedBed()
    req(selectedBed)
    data <- fread(selectedBed$path, sep = "\t", header = TRUE)
    if(ncol(data) == 3) {
        setnames(data, bed3Cols)
        data[, name := "."]
    } else if(ncol(data) == 4) {
        setnames(data, bed4Cols)
    } else {
        data <- data[, 1:5]
        setnames(data, c(bed4Cols, names(data)[5]))
    }
    setkeyv(data, bed3Cols)
    data
})

#----------------------------------------------------------------------
# summary metadata on the selected regions relative to the genome
#----------------------------------------------------------------------
regionsSummary <- reactive({
    selectedBedData <- selectedBedData()
    req(selectedBedData)
    sourceId <- sourceId()
    req(sourceId)
    ref <- getSourcePackageOption(sourceId, "composite", "primary-genome")
    gm <- genome_metadata[[ref]]
    regions_bp <- selectedBedData[, sum(end1 - start0)]
    list(
        total_bp           = gm$total_bp,
        included_bp        = gm$included_bp,
        regions_bp         = regions_bp,
        fraction_of_genome = regions_bp / gm$included_bp
    )
})
output$regionsSummary <- renderUI({
    summary <- regionsSummary()
    req(summary)
    tagList(
        tags$p("Total genome bp: ",    format(summary$total_bp,    big.mark = ",", scientific = FALSE)),
        tags$p("Included genome bp: ", format(summary$included_bp, big.mark = ",", scientific = FALSE)),
        tags$p("Total regions bp: ",   format(summary$regions_bp,  big.mark = ",", scientific = FALSE)),
        tags$p("Percent of included genome: ", format(summary$fraction_of_genome * 100, digits = 3), "%")
    )
})

#----------------------------------------------------------------------
# return reactives for data and metadata recovery
#----------------------------------------------------------------------
list(
    path = reactive({ selectedBed()$path }),
    name = reactive({ selectedBed()$name }),
    data = selectedBedData,
    metadata = regionsSummary
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
