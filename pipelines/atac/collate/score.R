# action:
#     calculate and analyze the distributions of the following score metrics per bin:
#         genome-level (not specific to samples at different spermiogenic stages):
#               gc   = aggregated bin GC content over all samples
#               txn  = (nascent) Transcription count per million reads
#         sample-level (specific to each sample at different spermiogenic stages):
#               gcrz = GC Residual Z-Score (excess or deficit of reads relative to GC peers)
#               nrll = protamine- vs. histone-associated Normalized Relative Log Likelihood
#     scores are always calculated relative to the primary genome distribution
#         scores may be present for the spike-in genome, but should generally be disregarded as unmeaningful
# input:
#     sample metadata file
#     composite genome bins BED file
#     insert_spans file created by atac/sites
#     output of the atac/collate/collate action step
# outputs:
#     a nested list with information all all scores above aggregated by:
#         individual samples
#         spermiogenic stages
#         difference between early and late spermiogenic stages
#     where each level has:
#         raw or aggregated scores per bin
#         the distribution of bin score
#         the mean and standard deviation of bin scores
#         the z-score of bin scores
#         the quantile of bin scores
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
        'METADATA_FILE',
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'GENOME_FASTA_SHM',
        'GENOME_BINS_BED',
        'INSERT_SPANS_DIR',
        'MAPPABILITY_SIZE_LEVELS',
        'ACTION_DIR',
        'DATA_FILE_PREFIX',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE',
        'STAGE_TYPES',
        'TRANSCRIPTION_BED',
        'SHM_FILE_PREFIX'
    ),
    integer = c(
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'BIN_SIZE',
        'N_CPU'
    )
))
if(env$TRANSCRIPTION_BED == "NA") env$TRANSCRIPTION_BED <- paste(env$DATA_FILE_PREFIX, "nascent_transcriptome_unstranded.bed.gz", sep = '.')
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$MODULES_DIR, 'bin')
sourceScripts(rUtilDir, c('bin_functions'))
rUtilDir <- file.path(env$MODULES_DIR, 'score')
sourceScripts(rUtilDir, c('score_functions'))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------
message("loading collate step output")
collate  <- readRDS(paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.'))

message("parsing spermiogenic stage types and GC limits")
stageTypes <- unpackStageTypes(env)
gcLimits <- unpackGcLimits(env)

message("calculating aggregated bin GC content")
collate$bins[, gc := rowMeans(collate$gc_bin_smp, na.rm = TRUE)]

message("extracting histone- and protamine-associated insert size distributions for primary genome")
emissProbsFile <- extractInsertSizeEps(collate$samples, collate$f_obs_ref_isl_smp$primary, env)

message("analyzing genome-level scores")
scores <- list(genome = list())
message("  gc")
scores$genome$gc <- analyzeScoreDist(
    collate$bins, 
    env$PRIMARY_GENOME, 
    gcLimits, 
    collate$bins$gc, 
    scoreTypes$genome$gc
)
scores$genome$gc$score <- NULL 
scores$genome$txn <- if(file.exists(env$TRANSCRIPTION_BED)) {
    message(paste("  txn:", env$TRANSCRIPTION_BED))
    txn_cpm <- fread(env$TRANSCRIPTION_BED)[[4]] # requires BED4 with normalized bin transcription value in column 4
    analyzeScoreDist(
        collate$bins, 
        env$PRIMARY_GENOME, 
        gcLimits, 
        txn_cpm, 
        scoreTypes$genome$txn
    )
} else {
    message("  no txn data")
    NULL
}

message("analyzing sample-level scores")
scores$sample <- sapply(names(scoreTypes$sample), function(scoreTypeName) {
    scoreType <- scoreTypes$sample[[scoreTypeName]]
    if(scoreType$gcBiasDependent) return(NULL)
    message(paste(" ", scoreTypeName))
    scoreFnName <- paste("get", scoreTypeName, sep = "_")
    sampleScores  <- analyzeSampleScores(
        collate,
        env$PRIMARY_GENOME, 
        gcLimits, 
        get(scoreFnName), 
        env, 
        scoreType, 
        emissProbsFile
    )
    aggregateScores <- aggregateSampleScores(
        collate,
        env$PRIMARY_GENOME, 
        gcLimits, 
        stageTypes, 
        sampleScores, 
        env, 
        scoreType
    )
    if(!("score" %in% scoreType$include)){
        for(key in names(sampleScores)) sampleScores[[key]]$score <- NULL
        for(key in names(aggregateScores$by_stage)) aggregateScores$by_stage[[key]]$score <- NULL
        for(key in names(aggregateScores$by_stageType)) aggregateScores$by_stageType[[key]]$score <- NULL
        aggregateScores$stageType_delta$score <- NULL 
    }
    list(
        sampleScores    = sampleScores,
        aggregateScores = aggregateScores
    )
}, simplify = FALSE, USE.NAMES = TRUE)

message()
message("saving output for app")
env$MAPPABILITY_SIZE_LEVELS <- as.integer(strsplit(env$MAPPABILITY_SIZE_LEVELS, " ")[[1]])
obj <- list(
    env = env[c(
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'MAPPABILITY_SIZE_LEVELS',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'BIN_SIZE',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE'
    )],
    samples     = collate$samples,
    references  = collate$references,
    stageTypes  = stageTypes,
    gcLimits    = gcLimits,
    scores      = scores
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "scores.rds", sep = '.')
)
str(obj)
#=====================================================================================
