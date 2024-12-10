# action:
#     in parallel over multiple samples:
#         assess the likelihood of bin insert sizes relative to two training states
#         segment the genome over all timed samples into nucleosomal and subnucleosomal output states
# input:
#     sample metadata file
#     genome and spike-in bins BED files
#     aligned, sorted, and indexed bam files (not necessarily deduplicated, we will dedup)
#     output of the atac/collate action
# outputs:
#     with data for each bin in each sample:
#        NRLL object = normalized relative log likelihoods = (LL_subnucleosomal - LL_nucleosomal) / nInserts
#        HMM object  = hidden Markov model (HMM) nucleosomal/subnucleosomal state output

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
        'SHM_FILE_PREFIX'
    ),
    integer = c(
        'MIN_MAPQ',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'BIN_SIZE',
        'N_CPU'
    ),
    double = c(
        'TRANSITION_PROBABILITY'
    )
))
persistence <- log(1 - env$TRANSITION_PROBABILITY)
# #-------------------------------------------------------------------------------------
# # source R scripts
# rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
# sourceScripts(rUtilDir, c('utilities','simple_HMM'))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------

message("loading collate action outputs")
# bd  <- readRDS(paste(env$DATA_FILE_PREFIX, "binCounts.rds",   sep = '.'))
isd <- readRDS(paste(env$DATA_FILE_PREFIX, "insertSizes.rds", sep = '.'))
nSamples <- nrow(isd$samples)

message("extracting nucleosomal and subnucleosomal insert size distributions")
getStateEmissProbs <- function(stage_){ # where a specific stage is taken as being a sufficiently pure representation of a state
    stage_samples <- isd$samples[stage == stage_, sample_name]
    x <- rowSums(isd$insertSizes$genome[, .SD, .SDcols = stage_samples])
    x <- x / sum(x)   # express as a proportion of the total
    x <- pmax(1e-5, x) # prevent log(0) and impossible values
    log(x / sum(x))   # normalize to sum to 1 and take the log for NRLL calculation and HMM
}
eps <- data.table(
    nucleosomal    = getStateEmissProbs(env$NUCLEOSOMAL_STAGE),
    subnucleosomal = getStateEmissProbs(env$SUBNUCLEOSOMAL_STAGE)
)
emissProbsFile <- paste(env$SHM_FILE_PREFIX, "emissionProbs.tsv", sep = '.')
write.table(
    eps, 
    file = emissProbsFile, 
    quote = FALSE, 
    row.names = FALSE, 
    col.names = FALSE, 
    sep = "\t"
)

message("calculating bin log likelihoods and solving HMM for subnucleosomal and nucleosomal states")
ref <- isd$references$genome
solveBinHMM <- paste('bash', file.path(env$ACTION_DIR, 'solve_bin_HMM.sh'))
solution <- do.call(cbind, mclapply(1:nSamples, function(sampleI) {
# solution <- do.call(cbind, lapply(1:nSamples, function(sampleI) {
    sample <- isd$samples[sampleI]
    message(paste('   ', sample$filename_prefix, '=', sample$sample_name))
    bamFile <- file.path(ref$input_dir, paste0(sample$filename_prefix, '.*.bam'))
    do.call(rbind, lapply(ref$chroms, function(chrom) { # all bins values over all ordered chroms
        fread(cmd = paste(solveBinHMM, bamFile, chrom, ref$fai_file, emissProbsFile, persistence)) # one value per bin on chrom
        # print(head(x[3000:100000], 100))
        # stop("stop")
    }))
}, mc.cores = env$N_CPU))
# }))

str(solution)

message("parsing HMM output states")
HMM <- do.call(cbind, mclapply(1:nSamples, function(sampleI) {
# HMM <- do.call(cbind, lapply(1:nSamples, function(sampleI) {
    stateI <- 1L + (sampleI - 1L) * 4L # solution dt has 4 columns per sample: state, LL_nuc, LL_subnuc, and nInserts
    data.table(
        solution[[stateI]] # encoded as 0 for nucleosomal and 1 for subnucleosomal states
    )
}, mc.cores = env$N_CPU))
# }))
names(HMM) <- isd$samples$sample_name

str(HMM)

message("parsing normalized relative log likelihoods")
NRLL <- do.call(cbind, mclapply(1:nSamples, function(sampleI) {
# NRLL <- do.call(cbind, lapply(1:nSamples, function(sampleI) {
    LL_nucI    <- 2L + (sampleI - 1L) * 4L # solution dt has 4 columns per sample: state, LL_nuc, LL_subnuc, and nInserts
    LL_subnucI <- LL_nucI + 1L
    nInsertsI  <- LL_nucI + 2L
    data.table(
        ifelse(
            solution[[nInsertsI]] > 0,
            (solution[[LL_subnucI]] - solution[[LL_nucI]]) / solution[[nInsertsI]],
            0
        )
    )
}, mc.cores = env$N_CPU))
# }))
names(NRLL) <- isd$samples$sample_name

str(NRLL)

# TODO: use rle to find HMM segment transitions and calculate NRLL for each alternating nuc/subnuc segment

message()
message("saving output for app")
obj <- list(
    bin_size    = env$BIN_SIZE,
    samples     = isd$samples,
    references  = isd$references,
    bins        = isd$bins,
    NRLL        = NRLL,
    HMM         = HMM
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "segmentation.rds", sep = '.')
)
str(obj)
#=====================================================================================
