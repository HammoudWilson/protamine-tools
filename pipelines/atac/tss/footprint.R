# action:
#     in parallel over multiple samples:
#         profile ATAC inserts at base-level resolution for the primary genome
#         create tabix-indexed files for footprint visualization
# input:
#     sample metadata file
#     insert_spans file created by atac/sites
# outputs:
#     tabix-indexed bgz files for every sample, stage, and stage type
#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'INSERT_SPANS_DIR',
        'MAPPABILITY_SIZE_LEVELS',
        'TN5_KMERS_FILE',
        'GENOME_METADATA_PREFIX',
        'ACTION_DIR',
        'DATA_FILE_PREFIX',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE',
        'STAGE_TYPES'
    ),
    integer = c(
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'N_CPU'
    ),
    double = c(
        'TOTAL_RAM_INT' # numeric since can exceed 2^31-1
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$MODULES_DIR, 'score')
sourceScripts(rUtilDir, c('score_functions'))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
getRamPerSort <- function(nSorts){
    floor((env$TOTAL_RAM_INT - 4e9) / nSorts)
}
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------
message("loading sample metadata")
samples <- fread(env$METADATA_FILE)[order(staging_order)]
# samples <- samples[filename_prefix %in% c("24290X11", "24290X9")]
# samples <- samples[filename_prefix %in% c("24290X8", "24290X10")]
nSamples <- nrow(samples)

message("parsing stages")
stages <- unique(samples$stage)
nStages <- length(stages)

message("parsing stage types")
stageTypes <- unpackStageTypes(env)
reverseStageTypes <- {
    reversed <- list()
    for (stageType in names(stageTypes)) {
        for (stage in stageTypes[[stageType]]) {
            reversed[[stage]] <- stageType
        }
    }
    reversed
}
nStageTypes <- length(stageTypes)

message("initializing filename object")
footprint <- list()

message("creating footprinting file by sample")
createFootprintFile <- file.path(env$ACTION_DIR, 'footprint_sample.sh')
bytesRamPerSort <- getRamPerSort(nSamples)
footprint$sample <- mclapply(1:nSamples, function(sampleI) {
# footprint$sample <- lapply(1:nSamples, function(sampleI) {
    smp <- samples[sampleI] 
    message(paste0('   ', smp$filename_prefix, ' = ', smp$sample_name, " (", smp$staging, ")"))
    stageType <- reverseStageTypes[[smp$stage]]
    if(is.null(stageType)) stageType <- 'NA'
    fread(cmd = paste("bash ", createFootprintFile, stageType, smp$stage, smp$sample_name, smp$filename_prefix, bytesRamPerSort))
}, mc.cores = env$N_CPU)
# })
names(footprint$sample) <- samples$sample_name

message("creating footprinting file by stage")
createFootprintFile <- file.path(env$ACTION_DIR, 'footprint_stage.sh')
bytesRamPerSort <- getRamPerSort(nStages)
footprint$stage <- mclapply(stages, function(stage) {
# footprint$stage <- lapply(stages, function(stage) {
    message(paste0('   ', stage))
    stageType <- reverseStageTypes[[stage]]
    if(is.null(stageType)) stageType <- 'NA'
    fread(cmd = paste("bash ", createFootprintFile, stageType, stage, bytesRamPerSort))
}, mc.cores = env$N_CPU)
# })
names(footprint$stage) <- stages

message("creating footprinting file by stage type")
createFootprintFile <- file.path(env$ACTION_DIR, 'footprint_stage_type.sh')
bytesRamPerSort <- getRamPerSort(nStageTypes)
footprint$stageType <- mclapply(names(stageTypes), function(stageType) {
# footprint$stageType <- lapply(names(stageTypes), function(stageType) {
    message(paste0('   ', stageType))
    fread(cmd = paste("bash ", createFootprintFile, stageType, bytesRamPerSort))
}, mc.cores = env$N_CPU)
# })
names(footprint$stageType) <- names(stageTypes)

message()
message("saving footprint file metadata for app")
env$INSERTS_BGZ_DIR <- file.path(env$TASK_DIR, "inserts_bgz")
obj <- list(
    env = env[c(
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE',
        'INSERTS_BGZ_DIR'
    )],
    samples             = samples,
    stages              = stages,
    stageTypes          = stageTypes,
    reverseStageTypes   = reverseStageTypes,
    footprint           = footprint
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "footprint.rds", sep = '.')
)
str(obj)
#=====================================================================================
