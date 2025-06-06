# action:
#     in parallel over multiple samples:
#         establish insert size distributions for each genome stratified by insert percent GC
#         count reads in composite genome bins, including primary and spike-in genomes
#         account for Tn5 cleavage site preference by calculating counting weights
#         calculate sample-specific bin mappability and GC content based on insert size distributions
# input:
#     sample metadata file
#     composite genome bins BED file
#     insert_spans file created by atac/sites
# outputs:
#     see below
# variable component abbreviations:
#     ref = a reference genome, i.e., primary, spike-in, or composite genome
#     smp = a sample, i.e., values representing an aggregated libary
#     ins = insert, i.e., a productively aligned read pair from a sample
#     is  = insert size, i.e., the reference genome span of an aligned read pair
#     isl = insert size level, i.e., the base insert size representing a range of insert sizes, e.g., isl = 35 for is[30, 35)
#     n   = number of inserts, i.e., a count of read pair observations
#     f   = a fraction of a total, i.e., a frequency
#     obs = observed n, or another value, prior to any normalization or negative binomial regression
#     exp = expected n, i.e., the predicted bin count after normalization and negative binomial regression
#     wgt = a site weight, used to adjust for Tn5 cleavage site preferences
#     chr = a chromosome, i.e., a reference genome sequence
#     bin = a 1kb genome bin on a chromosome
#     bs0 = bin start0, i.e., 0-based bin start position on a chromosome; paired chr and bs0 values uniquely identify a bin
#     gc  = fraction GC of a span of bases, either an insert or a bin
#     mpp = mappability of a span of bases, e.g., a bin
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
        'GENOME_METADATA_PREFIX',
        'TMP_FILE_PREFIX',
        'TN5_KMERS_FILE'
    ),
    integer = c(
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
wgtTypes <- c("observed", "weighted")

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
#  $ gc_35   : num  0 0 0 0 0 0 0 0 0 0 ... etc.
j <- seq(5, ncol(bins), 2)
mpp_isl_bin <- t(as.matrix(bins[, ..j])) # isl in rows, bins in columns, to apply matrix multiplication * f_smp_isl
j <- seq(6, ncol(bins), 2)
gc_isl_bin <- t(as.matrix(bins[, ..j]))

message("updating references based on bins")
for(refType in refTypes) {
    genome <- refs[[refType]]$genome
    I <- bins[, endsWith(chr, genome)]
    refs[[refType]]$chroms <- bins[I, unique(chr)]
    refs[[refType]]$n_bins <- sum(I)
    refs[[refType]]$binI <- range(which(I)) 
}

message("collating (merged) read pair data by sample")
collateSampleInserts <- paste('bash', file.path(env$ACTION_DIR, 'collate_sample_inserts.sh'))
sampleData <- mclapply(1:nrow(samples), function(sampleI) {
# sampleData <- lapply(1:nrow(samples), function(sampleI) {
    smp <- samples[sampleI]
    label <- paste0('   ', smp$filename_prefix, ' = ', smp$sample_name, " (", smp$staging, ")")
    message(label)

    # collate sample, recovering usable individual inserts with bin, insert size level, and site weight per insert
    # table has inserts for both primary and spike-in genomes, including sex chromosomes
    # ins_bin_isl sort is not guaranteed, but not required
    ins_bin_isl <- fread(cmd = paste(collateSampleInserts, smp$filename_prefix), sep = "\t", header = FALSE)
    setnames(ins_bin_isl, c('chr', 'bs0', 'isl','wgt'))

    # calculate the distribution of insert Tn5 weights (do this pre-limit to support plotting in app)
    # value is passed as int(round(log10(wgt) * 10))
    n_ins_wgt <- ins_bin_isl[wgt > 0, .(log10_wgt = as.integer(round(log10(wgt) * 10)))][, .(n_ins = .N), keyby = .(log10_wgt)]

    # enforce limits on weights to avoid excessive impact of extreme values, especially at high weights from low-preference sites
    # empirical observation indicates the following approximate quantiles across all samples:
    #     1%  ~ -2.2, i.e., 0.0063, >100-fold lower  weight than the unit insert count
    #     50% ~ -0.8, i.e., 0.16,   ~6-fold   lower  weight than the unit insert count (remember, there are a lot of these inserts)
    #     99% ~  1.1, i.e., 12.6,   ~12-fold  higher weight than the unit insert count (rare inserts have pre-limit values approaching 1e5!)
    #     ~10% of inserts have weights > 1
    min_wgt <- quantile(ins_bin_isl$wgt, 0.01, na.rm = TRUE)
    max_wgt <- quantile(ins_bin_isl$wgt, 0.99, na.rm = TRUE)
    ins_bin_isl[wgt < min_wgt, wgt := min_wgt]
    ins_bin_isl[wgt > max_wgt, wgt := max_wgt]

    # recover the insert size vs. gc distributions for the sample, stratified by genome and weighting
    # used for plotting in app
    n_is_gc <- fread(paste(env$TMP_FILE_PREFIX, smp$filename_prefix, "n_is_gc.txt", sep = "."), sep = "\t", header = FALSE)
    refs_is_gc <- n_is_gc[[1]]
    wgts_is_gc <- n_is_gc[[2]]
    n_is_gc <- n_is_gc[, -1:-2] # four sets of insert size vs. GC counts, one set per genome+weighting
    n_ref_wgt_is_gc <- sapply(refTypes, function(refType) {
        sapply(wgtTypes, function(weighting) {
            x <- as.matrix( n_is_gc[
                refs_is_gc == refs[[refType]]$genome &
                wgts_is_gc == weighting
            ])
            colnames(x) <- NULL
            x
        }, simplify = FALSE, USE.NAMES = TRUE)
    }, simplify = FALSE, USE.NAMES = TRUE)

    # count inserts per bin to establish observed bin counts
    # sum insert weights per bin to establish weighted bin counts
    # wgt values are properly calculated by collate_sample_inserts.pl using genome-specific tn5_site_freqs_exp
    bins_smp <- merge(
        bins[, .(chr, bs0)],
        ins_bin_isl[, .(n_obs = .N, n_wgt = sum(wgt)), keyby = .(chr, bs0)],
        all.x = TRUE,
        sort = FALSE,
        by = c('chr', 'bs0')
    )
    bins_smp[is.na(n_obs), ":="(n_obs = 0L, n_wgt = 0)]

    # count inserts by insert size level to establish weights for calculating sample-specific bin mappability and GC content
    # this is done using observed, not weighted, insert counts, since amplification and alignment bias were applied to 
    # the DNA fragments Tn5 actually gave us, not the theoretical reads that an unbiased transposase would have given us
    # it is also done per genome, since primary and spike-in genomes may have different insert size distributions
    f_obs_ref_isl <- sapply(refTypes, function(refType) {
        genome <- refs[[refType]]$genome
        x <- merge( # ensure that all insert size levels are present
            data.table(isl = as.integer(strsplit(env$MAPPABILITY_SIZE_LEVELS, " ")[[1]])),
            ins_bin_isl[endsWith(chr, genome), .(n_obs = .N), keyby = .(isl)],
            by = 'isl',
            all.x = TRUE,
            sort = FALSE
        )
        x[is.na(n_obs), n_obs := 0L]
        x[, n_obs / sum(n_obs)]
    }, simplify = FALSE, USE.NAMES = TRUE)

    # use the sample insert size level weights to calculate sample-specific bin mappability and GC content
    for(refType in refTypes) {
        genome <- refs[[refType]]$genome
        I <- bins[, endsWith(chr, genome)]
        bins_smp[I, ":="(
            mpp = colSums(mpp_isl_bin[, I] * f_obs_ref_isl[[refType]]),
            gc  = colSums(gc_isl_bin[, I]  * f_obs_ref_isl[[refType]])
        )]
    }

    # return sample-level data
    list(
        bins = bins_smp[, .(n_obs, n_wgt, mpp, gc)],
        f_obs_ref_isl   = f_obs_ref_isl,
        n_ref_wgt_is_gc = n_ref_wgt_is_gc,
        n_ins_wgt       = n_ins_wgt
    )
}, mc.cores = env$N_CPU)
# })
names(sampleData) <- samples$sample_name

message("refactoring sample data types into matrices (row=bin x col=sample)")
n_obs_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, n_obs]
})
colnames(n_obs_bin_smp) <- samples$sample_name
n_wgt_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, n_wgt]
})
colnames(n_wgt_bin_smp) <- samples$sample_name
mpp_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, mpp]
})
colnames(mpp_bin_smp) <- samples$sample_name
gc_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, gc]
})
colnames(gc_bin_smp) <- samples$sample_name
f_obs_ref_isl_smp <- sapply(refTypes, function(refType) {
    m <- sapply(samples$sample_name, function(sample_name) {
        sampleData[[sample_name]]$f_obs_ref_isl[[refType]]
    })
    colnames(m) <- samples$sample_name
    m
}, simplify = FALSE, USE.NAMES = TRUE)
n_ref_wgt_is_gc_smp <- sapply(refTypes, function(refType) {
    sapply(wgtTypes, function(weighting) {
        m1 <- sampleData[[1]]$n_ref_wgt_is_gc[[refType]][[weighting]]
        array(
            do.call(c, lapply(samples$sample_name, function(sample_name) {
                as.vector(sampleData[[sample_name]]$n_ref_wgt_is_gc[[refType]][[weighting]])
            })),
            dim = c(nrow(m1), ncol(m1), nSamples),
            dimnames = list(
                insert_size = NULL,
                gc_bin      = NULL,
                sample      = samples$sample_name
            )
        )
    }, simplify = FALSE, USE.NAMES = TRUE)
}, simplify = FALSE, USE.NAMES = TRUE)
n_ins_wgt_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$n_ins_wgt
}, simplify = FALSE, USE.NAMES = TRUE)

message()
message("saving collated sample data for app")
rm(sampleData)
env$MAPPABILITY_SIZE_LEVELS <- as.integer(strsplit(env$MAPPABILITY_SIZE_LEVELS, " ")[[1]])
obj <- list(
    env = env[c(
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'MAPPABILITY_SIZE_LEVELS',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'BIN_SIZE'
    )],
    samples    = samples,
    references = refs,
    bins = bins[, .(
        chrom = chr, start0 = bs0, 
        included = included
    )],
    n_obs_bin_smp = n_obs_bin_smp,
    n_wgt_bin_smp = n_wgt_bin_smp,
    mpp_bin_smp   = mpp_bin_smp,
    gc_bin_smp    = gc_bin_smp,
    f_obs_ref_isl_smp   = f_obs_ref_isl_smp,
    n_ref_wgt_is_gc_smp = n_ref_wgt_is_gc_smp,
    n_ins_wgt_smp = n_ins_wgt_smp
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.')
)
str(obj)
