#----------------------------------------------------------------------
# data recovery from cuttag data packages
#----------------------------------------------------------------------
paCutTag_ttl <- CONSTANTS$ttl$month
paCutTag_force <- FALSE
paCutTag_create <- "asNeeded"
paCutTag_loadPersistent <- function(..., sep = "\t", header = FALSE, force = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(force)) force <- paCutTag_force
    filePath <- loadPersistentFile(..., sep = sep, header = header, force = force, ttl = paCutTag_ttl)
    persistentCache[[filePath]]$data
}
paCutTag_getCached <- function(..., create = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(create)) create <- paCutTag_create
    protaminerCache$get(..., permanent = TRUE, create = create)$value
}
paCutTag_break_file <- function(sourceId, type){
    expandSourceFilePath(sourceId, paste0("cuttag-", type, ".rds"))
} 
paCutTag_get_component <- function(sourceId, type){
    typeFile <- paCutTag_break_file(sourceId, type)
    if(file.exists(typeFile)){
        readRDS(typeFile)
    } else {
        startSpinner(session, message = "loading cuttag data")
        x <- readRDS(getSourceFilePath(sourceId, "cuttag"))
        for(type_ in names(x)){
            startSpinner(session, message = paste("breaking", type_))
            filePath <- paCutTag_break_file(sourceId, type_)
            if(!file.exists(filePath)) saveRDS(x[[type_]], file = filePath)
        }
        x[[type]]
    }
}
paCutTag_load_ram <- function(sourceId, type){
    paCutTag_getCached(
        paste("cuttag", type, sep = "-"),
        key = sourceId,
        from = 'ram',
        createFn = function(...) {
            paCutTag_get_component(sourceId, type)
        },
        spinnerMessage = paste("loading", type)
    )
}

# cuttag samples sorted by staging_order
paCutTag_samples <- function(sourceId) paCutTag_load_ram(sourceId, "samples")[order(antibody_target, staging_order)]
paCutTag_samples_reactive <- function(sourceId) reactive({
    sourceId <- sourceId()
    req(sourceId)
    paCutTag_samples(sourceId)
})

#----------------------------------------------------------------------
# score metadata retrieval (not the actual scores)
#----------------------------------------------------------------------

# load score summaries and distributions into RAM
paCutTag_metadata <- function(sourceId){
    paCutTag_getCached(
        "cuttag_summary",
        key = sourceId,
        from = 'ram',
        createFn = function(...) {
            readRDS(getSourceFilePath(sourceId, "cuttag"))
        },
        spinnerMessage = paste("loading cuttag summary")
    )
}
getCutTagSampleMetadataList <- function(sourceId, scoreTypeName){ # returns a list of sample-level score objects based on GC normalization
    paCutTag_metadata(sourceId)$scores$sample[[scoreTypeName]]
}
getCutTagSampleMetadata <- function(sourceId, scoreTypeName, samples){ # returns a list of sample-level score objects
    getCutTagSampleMetadataList(sourceId, scoreTypeName)$sampleScores[samples$sample_name]
}
getCutTagStageMetadata <- function(sourceId, scoreTypeName, samples){ # returns a list of stage-level score objects matching a list of samples
    getCutTagSampleMetadataList(sourceId, scoreTypeName)$aggregateScores$by_stage[unique(samples$stage)]
}
