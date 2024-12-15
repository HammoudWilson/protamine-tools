# action:
#     calculate and analyze the distributions of the following score metrics per bin:
#         genome-level (not specific to samples at different spermiogenic stages):
#               gc   = bin GC content
#               txn  = (nascent) Transcription count per million reads
#         sample-level (specific to each sample at different spermiogenic stages):
#               gcrz = GC Residual Z-Score (excess or deficit of reads relative to GC peers)
#               iisf = Intermediate Insert Size Fraction (fraction of bin inserts with sizes >= 65 and <= 125bp)
#               nrll = protamine- vs. histone-associated Normalized Relative Log Likelihood
# input:
#     sample metadata file
#     genome and spike-in bins BED files
#     aligned, sorted, and indexed bam files (not necessarily deduplicated, we will dedup)
#     output of the atac/collate/collate action step
# outputs:
#     a nested list with information all all score above aggregated by:
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
        'GENOME_INPUT_DIR',
        'SPIKE_IN_INPUT_DIR',
        'GENOME',
        'SPIKE_IN_GENOME',
        'GENOME_FASTA',
        'SPIKE_IN_FASTA',
        'GENOME_BINS_BED',
        'SPIKE_IN_BINS_BED',
        'ACTION_DIR',
        'DATA_FILE_PREFIX',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE',
        'STAGE_TYPES',
        'TRANSCRIPTION_BED',
        'SHM_FILE_PREFIX'
    ),
    integer = c(
        'MIN_MAPQ',
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

message("loading collate step outputs")
bd  <- readRDS(paste(env$DATA_FILE_PREFIX, "binCounts.rds",   sep = '.'))
isd <- readRDS(paste(env$DATA_FILE_PREFIX, "insertSizes.rds", sep = '.'))
nSamples <- nrow(isd$samples)

message("parsing spermiogenic stage types and GC limits")
stageTypes <- unpackStageTypes(env)
gcLimits <- unpackGcLimits(env)

message("extracting histone- and protamine-associated insert size distributions")
emissProbsFile <- extractInsertSizeEps(isd, env)

message("analyzing genome-level scores")
scores <- list(genome = list())
message("  gc")
scores$genome$gc <- analyzeScoreDist(bd$bins$genome, gcLimits, bd$bins$genome$pct_gc, scoreTypes$genome$gc)
scores$genome$gc$score <- NULL 
scores$genome$txn <- if(file.exists(env$TRANSCRIPTION_BED)) {
    message(paste("  txn:", env$TRANSCRIPTION_BED))
    txn_cpm <- fread(env$TRANSCRIPTION_BED)[[4]] # requires BED4 with normalized bin transcription value in column 4
    analyzeScoreDist(bd$bins$genome, gcLimits, txn_cpm, scoreTypes$genome$txn)
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
    sampleScores  <- analyzeSampleScores(bd, gcLimits, get(scoreFnName), env, scoreType, emissProbsFile)
    aggregateScores <- aggregateSampleScores(bd, gcLimits, stageTypes, sampleScores, env, scoreType)
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
obj <- list(
    env         = env[c("BIN_SIZE","MAX_INSERT_SIZE","GENOME","HISTONE_STAGE","PROTAMINE_STAGE")],
    samples     = isd$samples,
    references  = isd$references,
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
