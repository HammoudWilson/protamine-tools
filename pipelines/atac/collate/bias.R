# action:
#     in parallel over multiple samples:
#         pull the 19-mer sequence surrounding the Tn5 cleavage site on forward alignment 5' ends
#         assemble a logo plot of the 19-mers per sample
# input:
#     sample metadata file
#     genome fa and fai files
#     aligned bam files
# outputs:
#     logo plot per sample per genome (grid of 19-mers by 4 bases)

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
    library(ggseqlogo)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'METADATA_FILE',
        'INPUT_DIR',
        'GENOME',
        'SPIKE_IN_GENOME',
        'GENOME_FASTA',
        'SPIKE_IN_FASTA',
        'GENOME_BINS_BED',
        'SPIKE_IN_BINS_BED',
        'GENOME_EXCLUSIONS_BED',
        'SPIKE_IN_EXCLUSIONS_BED',
        'ACTION_DIR',
        'DATA_FILE_PREFIX',
        'BIAS_PLOT_DIR'
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
# options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------
message("loading sample metadata")
samples <- fread(env$METADATA_FILE)
####################
# samples <- samples[filename_prefix %in% c("24290X4", "24290X5")] # one early_RS, one int_ES
nSamples <- nrow(samples)

message("parsing references")
references <- list(
    genome = list(
        folder_type = 'genome_folder',
        genome      = env$GENOME,
        fa_file     = env$GENOME_FASTA,
        bins_bed    = env$GENOME_BINS_BED,
        excl_bed    = env$GENOME_EXCLUSIONS_BED
    ),
    spike_in = list(
        folder_type = 'spike_in_folder',
        genome      = env$SPIKE_IN_GENOME,
        fa_file     = env$SPIKE_IN_FASTA,
        bins_bed    = env$SPIKE_IN_BINS_BED,
        excl_bed    = env$SPIKE_IN_EXCLUSIONS_BED
    )
)
refTypes <- names(references)

message("extracting Tn5 site 19-mers")
getLogoData <- paste('bash', file.path(env$ACTION_DIR, 'get_logo_data.sh'))
sapply(refTypes, function(refType) { # so, one list entry for genome, one for spike-in, same as bins
    message(paste(' ', refType))
    ref <- references[[refType]]
    mclapply(1:nSamples, function(sampleI) {
    # lapply(1:nSamples, function(sampleI) {
        sample <- samples[sampleI]
        sampleLabel <- paste0(sample$filename_prefix, ' ', sample$sample_name, ' (', sample$staging, ')')
        message(sampleLabel)
        bamFile <- file.path(env$INPUT_DIR, sample$base_folder, sample[[ref$folder_type]], paste0(sample$filename_prefix, '.*.bam'))
        
        # extract base frequencies for all 19-mer bases based on the 5' end of forward aligned reads
        logoData <- as.matrix(fread(
            cmd = paste(getLogoData, bamFile, ref$genome, ref$fa_file, ref$excl_bed, sample$filename_prefix, sample$staging_order),
            header = FALSE
        ))
        row.names(logoData) <- c('A', 'C', 'G', 'T')

        # create a logo plot of the 19-mer bases
        paddedStagingOrder <- sprintf("%02d", sample$staging_order)
        logoPlotFile <- file.path(env$BIAS_PLOT_DIR, paste0(paddedStagingOrder, '.', refType, '.', sample$filename_prefix, '.logo.png'))
        png(
            filename = logoPlotFile,
            width = 3.9,
            height = 1.9,
            units = 'in',
            res = 600
        )
        print( suppressWarnings(ggseqlogo(logoData)) )
        dev.off()

        # report the information content at the 6th position, the most informative position across all samples
        H6 <- -sum(sapply(1:4, function(i) logoData[i,6] * log2(logoData[i,6]))) # base uncertainty at the 6th position
        R6 <- log2(4) - H6 # information content at the 6th position
        message(paste(sample$staging_order, sample$staging, sample$stage, refType, sample$filename_prefix, sample$batch, "R6", R6, sep = '\t'))

    }, mc.cores = env$N_CPU)
    # })
})

# message()
# message("saving binCounts for app")
# obj <- list(
#     bin_size    = env$BIN_SIZE,
#     samples     = samples,
#     references  = references,
#     bins        = bins,
#     binCounts   = binCounts
# )
# saveRDS(
#     obj, 
#     file = paste(env$DATA_FILE_PREFIX, "binCounts.rds", sep = '.')
# )
# str(obj)

# message()
# message("saving insertSizes for app")
# obj <- list(
#     bin_size    = env$BIN_SIZE,
#     samples     = samples,
#     references  = references,
#     insertSizes = insertSizes
# )
# saveRDS(
#     obj, 
#     file = paste(env$DATA_FILE_PREFIX, "insertSizes.rds", sep = '.')
# )
# str(obj)
#=====================================================================================
