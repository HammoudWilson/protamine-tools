# action:
#     in parallel over multiple samples:
#         count reads in genome and spike-in bins
#         establish insert size distributions
# input:
#     sample metadata file
#     genome and spike-in bins BED files
#     aligned, sorted, and indexed bam files (not necessarily deduplicated, we will dedup)
# outputs:
#     samples     = data.table of sample metadata
#     bins        = data.table of genome bins, same bin order as binCounts
#     binCounts   = array of read counts in each bin for each sample, stratified by all_inserts and intermediate insert sizes
#     insertSizes = data.table of insert size distributions for each sample
#     where bins, binCounts, and insertSizes are named lists with separate entries for genome and spike-in reference types

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
        'DATA_FILE_PREFIX'
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
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------

message("loading sample metadata")
samples <- fread(env$METADATA_FILE)
# samples <- samples[filename_prefix %in% c("24290X11", "24290X9")]
nSamples <- nrow(samples)

message("parsing references")
references <- list(
    genome = list(
        input_dir   = env$GENOME_INPUT_DIR,
        genome      = env$GENOME,
        fai_file    = paste0(env$GENOME_FASTA, '.fai'),
        bins_bed    = env$GENOME_BINS_BED
    ),
    spike_in = list(
        input_dir   = env$SPIKE_IN_INPUT_DIR,
        genome      = env$SPIKE_IN_GENOME,
        fai_file    = paste0(env$SPIKE_IN_FASTA, '.fai'),
        bins_bed    = env$SPIKE_IN_BINS_BED
    )
)
refTypes <- names(references)

message("loading (spike-in) genome bins")
bins <- sapply(refTypes, function(refType) {
    bins <- fread(references[[refType]]$bins_bed) # already restricted by genome/bin to autosomes and chrX/Y
    references[[refType]]$chroms <<- bins[, unique(chrom)]
    references[[refType]]$n_bins <<- nrow(bins)
    bins[, nAlleles := ifelse(chrom %in% c('chrX', 'chrY'), 1L, 2L)]
    bins
}, simplify = FALSE, USE.NAMES = TRUE)

message("calculating sample (spike-in) bin counts")
countChromBins <- paste('bash', file.path(env$ACTION_DIR, 'count_chrom_bins.sh'))
insertTypes <- c('all_inserts','intermediate')
nInsertTypes <- length(insertTypes)
binCounts <- sapply(refTypes, function(refType) { # so, one list entry for genome, one for spike-in, same as bins
    message(paste(' ', refType))
    ref <- references[[refType]]
    array( # each refType is an array with dim1 = bins, dim2 = insertTypes, dim3 = samples
        do.call(c, mclapply(1:nSamples, function(sampleI) {
        # do.call(c, lapply(1:nSamples, function(sampleI) {
            sample <- samples[sampleI]
            message(paste('   ', sample$filename_prefix, '=', sample$sample_name))
            bamFile <- file.path(ref$input_dir, paste0(sample$filename_prefix, '.*.bam'))
            unlist(lapply(insertTypes, function(insertType) { # bins concatenated in two chunks, one per insert type
                unlist(lapply(ref$chroms, function(chrom) { # all bins values over all ordered chroms
                    fread(cmd = paste(countChromBins, bamFile, insertType, chrom, ref$fai_file)) # one value per bin on chrom
                }))
            }))
        }, mc.cores = env$N_CPU)),
        # })),
        dim = c(
            ref$n_bins, 
            nInsertTypes, 
            nSamples
        ), 
        dimnames = list(
            bin = NULL, 
            insert_type = insertTypes, 
            sample_name = samples$sample_name
        )
    )
}, simplify = FALSE, USE.NAMES = TRUE)

message("processing insert size distributions")
getInsertSizes <- paste('bash', file.path(env$ACTION_DIR, 'get_insert_sizes.sh'))
insertSizes <- sapply(refTypes, function(refType) { # so, one list entry for genome, one for spike-in, same as bins
    message(paste(' ', refType))
    ref <- references[[refType]]
    insertSizes <- do.call(cbind, mclapply(1:nrow(samples), function(sampleI) {
    # insertSizes <- do.call(cbind, lapply(1:nrow(samples), function(sampleI) {
        sample <- samples[sampleI]
        message(paste('   ', sample$filename_prefix, '=', sample$sample_name))
        bamFile <- file.path(ref$input_dir, paste0(sample$filename_prefix, '.*.bam'))
        fread(cmd = paste(getInsertSizes, bamFile))
    }, mc.cores = env$N_CPU))
    # }))
    setnames(insertSizes, samples$sample_name)
    insertSizes
}, simplify = FALSE, USE.NAMES = TRUE)

message()
message("saving binCounts for app")
obj <- list(
    bin_size    = env$BIN_SIZE,
    samples     = samples,
    references  = references,
    bins        = bins,
    binCounts   = binCounts
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "binCounts.rds", sep = '.')
)
str(obj)

message()
message("saving insertSizes for app")
obj <- list(
    bin_size    = env$BIN_SIZE,
    samples     = samples,
    references  = references,
    insertSizes = insertSizes
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "insertSizes.rds", sep = '.')
)
str(obj)
#=====================================================================================
