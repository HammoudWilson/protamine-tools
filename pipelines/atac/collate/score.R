# action:
#     calculate and analyze the distributions of the following score metrics per bin:
#         genome-level (not specific to samples at different spermiogenic stages):
#               gc   = bin GC content
#               rpkm = (nascent) transcription Reads Per Kilobase per Million reads
#         sample-level (specific to each sample at different spermiogenic stages):
#               cpm  = read Counts Per Million reads
#               gcrz = GC Residual Z-Score (excess or deficit of reads relative to GC peers)
#               snif = Subnucleosomal Insert Fraction (fraction of bin insert sizes < 146bp)
#               nrll = subnucleosomal vs. nucleosomal Normalized Relative Log Likelihood
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
        'NUCLEOSOMAL_STAGE',
        'SUBNUCLEOSOMAL_STAGE',
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

message("extracting nucleosomal and subnucleosomal insert size distributions")
emissProbsFile <- extractInsertSizeEps(isd, env)

message("analyzing genome-level scores")
scores <- list(genome = list())
message("  GC")
scores$genome$gc <- analyzeScoreDist(bd$bins$genome, gcLimits, bd$bins$genome$pct_gc, scoreTypes$genome$gc$distUnit)
scores$genome$rpkm <- if(env$TRANSCRIPTION_BED != "NA") {
    message("  RPKM")
    stop("rpkm analysis not implemented yet")
    analyzeScoreDist(bd$bins$genome, gcLimits, xxxx, scoreTypes$genome$rpkm$distUnit)
} else {
    message("  RPKM: no transcription data")
    NULL
}

message("analyzing sample-level scores")
scores$sample <- sapply(names(scoreTypes$sample), function(scoreTypeName) {
    scoreType <- scoreTypes$sample[[scoreTypeName]]
    if(scoreType$gbBiasDependent) return(NULL)
    message(paste(" ", scoreTypeName))
    scoreFnName <- paste("get", scoreTypeName, sep = "_")
    sampleScores  <- analyzeSampleScores(bd, gcLimits, get(scoreFnName), env, scoreType$distUnit, emissProbsFile)
    aggregateScores <- aggregateSampleScores(bd, gcLimits, stageTypes, sampleScores, env, scoreType$distUnit)
    list(
        sampleScores    = sampleScores,
        aggregateScores = aggregateScores
    )
}, simplify = FALSE, USE.NAMES = TRUE)

message()
message("saving output for app")
obj <- list(
    env         = env[c("BIN_SIZE","MAX_INSERT_SIZE","GENOME","NUCLEOSOMAL_STAGE","SUBNUCLEOSOMAL_STAGE")],
    samples     = isd$samples,
    references  = isd$references,
    scoreTypes  = scoreTypes,
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
