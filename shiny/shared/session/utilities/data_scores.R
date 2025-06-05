# data recover from atac/collate data packages
paScores_ttl <- CONSTANTS$ttl$month
paScores_force <- FALSE
paScores_create <- "asNeeded"
paScores_loadPersistent <- function(..., sep = "\t", header = FALSE, force = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(force)) force <- paScores_force
    filePath <- loadPersistentFile(..., sep = sep, header = header, force = force, ttl = paScores_ttl)
    persistentCache[[filePath]]$data
}
paScores_getCached <- function(..., create = NULL, spinnerMessage = NULL){
    if (!is.null(spinnerMessage)) startSpinner(session, message = spinnerMessage)
    if (is.null(create)) create <- paScores_create
    paScoresCache$get(..., permanent = TRUE, create = create)$value
}

# break a data package into its component object types and return the requested type object
paScores_break_file <- function(sourceId, type){
    expandSourceFilePath(sourceId, paste0("scores-", type, ".rds"))
}   
paScores_get_component <- function(sourceId, type){
    typeFile <- paScores_break_file(sourceId, type)
    if(file.exists(typeFile)){
        readRDS(typeFile)
    } else {
        startSpinner(session, message = "loading scores data")
        x <- readRDS(getSourceFilePath(sourceId, "scores"))
        for(type_ in names(x)){
            startSpinner(session, message = paste("breaking", type_))
            filePath <- paScores_break_file(sourceId, type_)
            if(!file.exists(filePath)) saveRDS(x[[type_]], file = filePath)
        }
        x[[type]]
    }
}


# # load and format genome bins and read counts (from atat/collate action)
# paScores <- function(sourceId){
#     startSpinner(session, message = "loading scores")
#     filePath <- loadPersistentFile(
#         sourceId = sourceId, 
#         contentFileType = "scores", 
#         ttl = CONSTANTS$ttl$month, 
#         postProcess = function(sd){
#             startSpinner(session, message = "scores post-processing")
#             sd$reverseStageTypes <- {
#                 reversed <- list()
#                 for (name in names(sd$stageTypes)) {
#                     for (value in sd$stageTypes[[name]]) {
#                         reversed[[value]] <- name
#                     }
#                 }
#                 reversed
#             }
#             sd
#         }
#     )
#     stopSpinner(session)
#     persistentCache[[filePath]]$data
# }


    # env = env[c(
    #     'PRIMARY_GENOME',
    #     'SPIKE_IN_GENOME',
    #     'GENOME',
    #     'MAPPABILITY_SIZE_LEVELS',
    #     'MIN_INSERT_SIZE',
    #     'MAX_INSERT_SIZE',
    #     'BIN_SIZE',
    #     'HISTONE_STAGE',
    #     'PROTAMINE_STAGE'
    # )],
    # samples     = collate$samples,
    # references  = collate$references,
    # stageTypes  = stageTypes,
    # gcLimits    = gcLimits,
    # scores      = scores