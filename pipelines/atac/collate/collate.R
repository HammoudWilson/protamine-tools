# action:
#     in parallel over multiple samples:
#         establish insert size distributions for each genome stratified by insert percent GC
#         count reads in composite genome bins, including primary and spike-in genomes
# input:
#     sample metadata file
#     composite genome bins BED file
#     aligned, sorted, and indexed bam files
#       one per sample
#       aligned to composite genome
#       not deduplicated, we will dedup
# outputs:
#     samples     = data.table of sample metadata
#     bins        = data.table of genome bins, same bin order as binCounts
#     binCounts   = array of read counts in each bin for each sample, stratified by all_inserts and intermediate insert sizes
#     insertSizes = data.table of insert size distributions for each sample
#     where bins, binCounts, and insertSizes are named lists with separate entries for genome and spike-in reference types
# variable abbreviations:
#     ins = insert, i.e., a productively aligned read pair
#     n   = number of inserts, i.e., a count of read pair observations
#     raw = observed n, or another value, prior to any normalization or negative binomial regression
#     exp = expected n, i.e., the predicted bin count after normalization or negative binomial regression
#     bin = a 1kb genome bin
#     chr = a chromosome, i.e., a reference genome sequence
#     bs0 = bin start0, i.e., 0-based bin start position on a chromosome; paired chr and bs0 values uniquely identify a bin
#     gc  = fraction GC of a span of bases
#     mpp = mappability of a span of bases, e.g., a bin
#     is  = insert size, i.e., the reference genome span of an aligned read pair
#     isl = insert size level, i.e., the base insert size representing a range of insert sizes, e.g., isl = 35 for is in [30, 35)
#     f   = a fraction of a total, i.e., a frequency
#     smp = a sample, i.e., values representing an aggregated libary
#     ref = a reference genome, i.e., primary, spike-in, or combined genome
# thus:
#     n_bin      = number of inserts in a bin
#     n_bin_isl  = number of inserts in a bin at a given insert size level
#     gc_bin_isl = fraction GC in a bin at a given insert size level
#     f_smp_isl  = fraction of inserts in a sample at a given insert size level
#     etc.
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
        'INPUT_DIR',
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'GENOME_FASTA_SHM',
        'GENOME_BINS_BED',
        'GENOME_INCLUSIONS_BED',
        'ACTION_DIR',
        'DATA_FILE_PREFIX',
        'TMP_FILE_PREFIX'
    ),
    integer = c(
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'MIN_MAPQ',
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
samples <- fread(env$METADATA_FILE)[order(staging_order)]
# samples <- samples[filename_prefix %in% c("24290X11", "24290X9")]
# samples <- samples[filename_prefix %in% c("24290X8", "24290X10")]
nSamples <- nrow(samples)

message("parsing references")
refs <- list(
    primary = list(
        genome = env$PRIMARY_GENOME
    ),
    spike_in = list(
        genome = env$SPIKE_IN_GENOME
    )
)
refTypes <- names(refs)

message("loading composite genome bins")
bins <- fread(env$GENOME_BINS_BED) # already restricted by genome/bin to autosomes and chrX/Y
binsNames <- names(bins)
binsNames[1:3] <- c('chr', 'bs0', 'be1') # for script naming consistency
setnames(bins, binsNames)
# Classes ‘data.table’ and 'data.frame':  2860977 obs. of  44 variables:
#  $ chr     : chr  "chr1-mm39" "chr1-mm39" "chr1-mm39" "chr1-mm39" ...
#  $ bs0     : int  0 1000 2000 3000 4000 5000 6000 7000 8000 9000 ...
#  $ be1     : int  1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 ...
#  $ included: int  0 0 0 0 0 0 0 0 0 0 ...
#  $ map_35  : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ gc_35   : num  0 0 0 0 0 0 0 0 0 0 ...
j <- seq(5, ncol(bins), 2)
mpp_bin_isl <- t(as.matrix(bins[, ..j])) # isl in rows, bins in columns, to apply matrix multiplication * f_smp_isl
j <- seq(6, ncol(bins), 2)
gc_bin_isl <- t(as.matrix(bins[, ..j]))

message("updating references based on bins")
for(refType in refTypes) {
    genome <- refs[[refType]]$genome
    I <- bins[, endsWith(chr, genome)]
    refs[[refType]]$chroms <- bins[I, unique(chr)]
    refs[[refType]]$n_bins <- sum(I)
    bins[I, nAlleles := ifelse(chr %in% paste0(c('chrX', 'chrY'), "-", genome), 1L, 2L)]
}

message("collating (merged) read pair data by sample")
collateSampleInserts <- paste('bash', file.path(env$ACTION_DIR, 'collate_sample_inserts.sh'))
sampleData <- mclapply(1:nrow(samples), function(sampleI) {
# sampleData <- lapply(1:nrow(samples), function(sampleI) {
    smp <- samples[sampleI]
    label <- paste0('   ', smp$filename_prefix, ' = ', smp$sample_name, " (", smp$staging, ")")
    message(label)
    bamFile <- file.path(env$INPUT_DIR, paste0(smp$filename_prefix, '.*.bam'))
    tmpPrefix <- paste0(env$TMP_FILE_PREFIX, ".collate.", smp$filename_prefix)

    # collate sample, recovering usable inserts with bin and insert size level metadata per insert
    ins_bin_isl <- fread(cmd = paste(collateSampleInserts, bamFile, tmpPrefix, smp$sample_name), sep = "\t", header = FALSE)
    setnames(ins_bin_isl, c('chr', 'bs0', 'isl'))

    # count inserts per bin to establish raw bin counts, the response variable for the negative binomial regression
    bins_smp <- merge(
        bins[, .(chr, bs0)],
        ins_bin_isl[, 
            .(n_raw = .N), 
            keyby = .(chr, bs0)
        ],
        all.x = TRUE,
        sort = FALSE,
        by = c('chr', 'bs0')
    )
    bins_smp[is.na(n_raw), n_raw := 0L]
    # message()
    # message("bins_smp")
    # str(bins_smp)

    # count inserts by insert size level to establish weights for calculating sample-specific bin mappability and GC content
    f_smp_isl <- ins_bin_isl[, 
        .(n = .N), 
        keyby = .(isl)
    ][, 
        n / sum(n)
    ]
    # message()
    # message("f_smp_isl")
    # print(f_smp_isl)

    # use the sample insert size level weights to calculate sample-specific bin mappability and GC content
    bins_smp[, ":="(
        mpp = colSums(mpp_bin_isl * f_smp_isl),
        gc  = colSums(gc_bin_isl  * f_smp_isl)
    )]
    # message()
    # message("bins_smp")
    # str(bins_smp)
    # str(bins_smp[mpp > 0])

    # recover the insert size vs. gc distribution for the sample, stratified by genome
    # used for plotting in app
    n_is_gc <- fread(paste0(tmpPrefix, ".n_is_gc.txt"), sep = "\t", header = FALSE)
    refs_is_gc <- n_is_gc[[1]]
    n_is_gc <- n_is_gc[, -1]
    n_is_gc_ref <- sapply(refTypes, function(refType) {
        n_is_gc <- as.matrix( n_is_gc[refs_is_gc == refs[[refType]]$genome] )
        colnames(n_is_gc) <- NULL
        n_is_gc
    }, simplify = FALSE, USE.NAMES = TRUE)
    # message()
    # message("n_is_gc_ref")
    # str(n_is_gc_ref)

    # recover the Tn5 kmers used by the sample to establish Tn5 preferences
    # TODO: this comes to us complete and pre-sorted, so we don't need to carry kmer in this object??
    # TODO: this can be done in perl...
    n_tn5_smp <- fread(paste0(tmpPrefix, ".n_tn5_smp.txt"), sep = "\t", header = FALSE)
    setnames(n_tn5_smp, c('kmer', 'n'))
    n_tn5_smp[, f := n / sum(n)]
    # message()
    # message("n_tn5_smp")
    # str(n_tn5_smp)
    # print(n_tn5_smp[order(f)])

    # TODO: make a call to collect isl-dependent Tn5 kmers and use tn5_smp to collect a weight Tn5 preference

    # TODO: use MASS::glm.nb to fit a negative binomial regression model to the bins_smp data

    # stop("XXXXXXXXXXXXXXXXXXXXXX")


    # collect collation counts by stage
    n_ins_stage <- fread(paste0(tmpPrefix, ".n_ins_stage.txt"))
    setnames(n_ins_stage, c('sample_name', 'stage', 'n'))

    # return sample-level data
    list(
        n_is_gc_ref = n_is_gc_ref,
        bins_smp    = bins_smp,
        f_smp_isl   = f_smp_isl,
        n_ins_stage = n_ins_stage
    )
}, mc.cores = env$N_CPU)
# })
names(sampleData) <- samples$sample_name

message()
message("saving collated sampleData for app")
obj <- list(
    bin_size    = env$BIN_SIZE,
    samples     = samples,
    references  = refs,
    bins        = bins[, .(chrom = chr, start0 = bs0, end1 = be1, included = included, nAlleles = nAlleles)],
    sampleData  = sampleData
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.')
)
str(obj)

message()
message("collation summary")
summary <- do.call(rbind, lapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$n_ins_stage
}))
print(dcast(summary, sample_name ~ stage, value.var = 'n', fun.aggregate = sum, fill = 0L))

# message("processing insert size distributions, stratified by genome and GC content")
# getInsertSizes <- paste('bash', file.path(env$ACTION_DIR, 'get_insert_sizes.sh'))
# insertSizes <- mclapply(1:nrow(samples), function(sampleI) {
# # insertSizes <- lapply(1:nrow(samples), function(sampleI) {
#     sample <- samples[sampleI]
#     message(paste0('   ', sample$filename_prefix, ' = ', sample$sample_name, " (", sample$staging, ")"))
#     bamFile <- file.path(env$INPUT_DIR, paste0(sample$filename_prefix, '.*.bam'))
#     dt <- fread(cmd = paste(getInsertSizes, bamFile), sep = "\t", header = FALSE)
#     refs <- dt[[1]]
#     dt <- dt[, -1]
#     sapply(refTypes, function(refType) {
#         ref <- references[[refType]]
#         m <- as.matrix( dt[refs == ref$genome] )
#         colnames(m) <- NULL
#         m
#     }, simplify = FALSE, USE.NAMES = TRUE)
# }, mc.cores = env$N_CPU)
# # })
# names(insertSizes) <- samples$sample_name

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
# message()
# str(obj)
# rm(obj)
# message()

# message("calculating sample-specific bin mappability based on insert size distribution")
# # rowSums insertSizes
# # insert size to levels
# # fraction per insert size level
# # weighted mean of level-specific bin mappability
# # use below to correct bin counts (divde by mappability, with min mappability threshold)

# message("calculating sample bin counts, stratified by genome")
# countChromBins <- paste('bash', file.path(env$ACTION_DIR, 'count_chrom_bins.sh'))
# insertTypes <- c('all_inserts','intermediate')
# nInsertTypes <- length(insertTypes)

# # best thing here will be to establish what the best output data structure is
# # not going to use intermediate insert assessment any longer (just use NRLL, they proved redundant)

# binCounts <- sapply(refTypes, function(refType) { # so, one list entry for genome, one for spike-in, same as bins
#     message(paste(' ', refType))
#     ref <- references[[refType]]
#     array( # each refType is an array with dim1 = bins, dim2 = insertTypes, dim3 = samples
#         do.call(c, mclapply(1:nSamples, function(sampleI) {
#         # do.call(c, lapply(1:nSamples, function(sampleI) {
#             sample <- samples[sampleI]
#             message(paste0('   ', sample$filename_prefix, ' = ', sample$sample_name, " (", sample$staging, ")"))
#             bamFile <- file.path(env$INPUT_DIR, sample$base_folder, sample[[ref$folder_type]], paste0(sample$filename_prefix, '.*.bam'))
#             unlist(lapply(insertTypes, function(insertType) { # bins concatenated in two chunks, one per insert type

#                 # this chroms isn't correct, can just run all chroms over both genome, just like bins map is
#                 # alternatively, need to split the bins object into two by genome

#                 unlist(lapply(ref$chroms, function(chrom) { # all bins values over all ordered chroms
#                     fread(cmd = paste(countChromBins, bamFile, insertType, chrom, ref$fai_file)) # one value per bin on chrom
#                 }))
#             }))
#         }, mc.cores = env$N_CPU)),
#         # })),
#         dim = c(
#             ref$n_bins, 
#             nInsertTypes, 
#             nSamples
#         ), 
#         dimnames = list(
#             bin = NULL, 
#             insert_type = insertTypes, 
#             sample_name = samples$sample_name
#         )
#     )
# }, simplify = FALSE, USE.NAMES = TRUE)

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
#=====================================================================================
