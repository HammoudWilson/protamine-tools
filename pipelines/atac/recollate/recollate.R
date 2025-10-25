# action:
#     in parallel over multiple samples, for primary genome only:
#         establish insert size distributions, stratified by insert percent GC and peak/non-peak
#         count reads in bins, now marking bins crossing peaks as not included
#         inherit sample-specific bin mappability and GC content from collate
# input:
#     sample metadata file
#     composite genome bins BED file
#     insert_spans file created by atac/sites
#     output of atac/collate
# outputs:
#     see below
# variable component abbreviations:
#     smp = a sample, i.e., values representing an aggregated libary
#     ins = insert, i.e., a productively aligned read pair from a sample
#     is  = insert size, i.e., the reference genome span of an aligned read pair
#     isl = insert size level, i.e., the base insert size representing a range of insert sizes, e.g., isl = 35 for is[30, 35)
#     n   = number of inserts, i.e., a count of read pair observations
#     f   = a fraction of a total, i.e., a frequency
#     exp = expected n, i.e., the predicted bin count after normalization and negative binomial regression
#     chr = a chromosome, i.e., a reference genome sequence
#     bin = a 1kb genome bin on a chromosome
#     bs0 = bin start0, i.e., 0-based bin start position on a chromosome; paired chr and bs0 values uniquely identify a bin
#     gc  = fraction GC of a span of bases, either an insert or a bin
#     mpp = mappability of a span of bases, e.g., a bin
#     peak= a called ab-initio accessibility peak
# unused variable component abbreviations (dropped relative to collate):
#     ref = a reference genome, i.e., primary, spike-in, or composite genome (recollate primary genome only)
#     obs = observed n, or another value, prior to any normalization or negative binomial regression (recollate observed counts only)
#     wgt = a site weight, used to adjust for Tn5 cleavage site preferences
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
        'TMP_FILE_PREFIX'
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

# message("loading collate output")
# collate <- readRDS(paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.'))
message("loading sample metadata")
samples <- fread(env$METADATA_FILE)[order(staging_order)]
nSamples <- nrow(samples)

message("loading composite genome bins")
bins <- fread(env$GENOME_BINS_BED) # already restricted by genome/bin to autosomes and chrX/Y

message("loading ab initio output")
abInitio <- readRDS(paste(env$DATA_FILE_PREFIX, "ab_initio.rds", sep = '.'))
regions <- abInitio$regions[index_stage == "overlap_group", .(chrom, start0, end1)]
regions[, chrom := paste0(chrom, '-', env$PRIMARY_GENOME)] 
rm(abInitio)
env$PEAK_REGIONS_BED <- paste(env$TMP_FILE_PREFIX, "peak_regions.bed", sep = '.')
fwrite(
    regions, 
    file = env$PEAK_REGIONS_BED, 
    col.names = FALSE,
    sep = "\t"
)

message("excluding bins that overlap accessibility peaks")
setkey(regions, chrom, start0, end1)
ov <- foverlaps(
    bins[, .(chrom, start0, end1)], 
    regions, 
    mult = "first",
    nomatch = NA,
    which = TRUE
)
bins[!is.na(ov), included := 0L]
rm(regions, ov)

message("finishing bin parsing")
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

message("recollating (merged) read pair data by sample")
recollateSampleInserts <- paste('bash', file.path(env$ACTION_DIR, 'recollate_sample_inserts.sh'))
sampleData <- mclapply(1:nrow(samples), function(sampleI) {
# sampleData <- lapply(1:nrow(samples), function(sampleI) {
    smp <- samples[sampleI]
    label <- paste0('   ', smp$filename_prefix, ' = ', smp$sample_name, " (", smp$staging, ")")
    message(label)

    # collate sample, recovering usable individual inserts with bin and insert size level per insert
    # table has inserts for primary genome only, including sex chromosomes
    # ins_bin_isl sort is not guaranteed, but not required
    ins_bin_isl <- fread(cmd = paste(recollateSampleInserts, smp$filename_prefix, env$PEAK_REGIONS_BED), sep = "\t", header = FALSE)
    setnames(ins_bin_isl, c('chr', 'bs0', 'isl', 'peak'))

    # recover the insert size vs. gc distributions for the sample, stratified by peak/non-peak
    # used for plotting in app
    n_is_gc <- fread(paste(env$TMP_FILE_PREFIX, smp$filename_prefix, "n_is_gc.txt", sep = "."), sep = "\t", header = FALSE)
    peak_is_gc <- n_is_gc[[1]]
    n_is_gc <- n_is_gc[, -1] # two sets of insert size vs. GC counts, one set per peak/non-peak
    n_peak_is_gc <- sapply(0:1, function(peak) {
        x <- as.matrix( n_is_gc[peak_is_gc == peak])
        colnames(x) <- NULL
        x
    }, simplify = FALSE, USE.NAMES = TRUE)

    # count inserts per bin to establish observed bin counts
    bins_smp <- merge(
        bins[, .(chr, bs0)],
        ins_bin_isl[, .(n_obs = .N), keyby = .(chr, bs0)],
        all.x = TRUE,
        sort = FALSE,
        by = c('chr', 'bs0')
    )
    bins_smp[is.na(n_obs), ":="(n_obs = 0L)]

    # count inserts by insert size level to establish weights for calculating sample-specific bin mappability and GC content
    # this is done using observed insert counts
    f_obs_isl_off_peak <- {
        x <- merge( # ensure that all insert size levels are present
            data.table(isl = as.integer(strsplit(env$MAPPABILITY_SIZE_LEVELS, " ")[[1]])),
            ins_bin_isl[peak == 0, .(n_obs = .N), keyby = .(isl)],
            by = 'isl',
            all.x = TRUE,
            sort = FALSE
        )
        x[is.na(n_obs), n_obs := 0L]
        x[, n_obs / sum(n_obs)]
    }

    # use the sample insert size level weights to calculate sample-specific bin mappability and GC content
    I <- bins[, endsWith(chr, env$PRIMARY_GENOME)]
    bins_smp[I, ":="(
        mpp = colSums(mpp_isl_bin[, I] * f_obs_isl_off_peak),
        gc  = colSums(gc_isl_bin[, I]  * f_obs_isl_off_peak)
    )]

    # return sample-level data
    list(
        bins = bins_smp[, .(n_obs, mpp, gc)],
        f_obs_isl    = f_obs_isl_off_peak,
        n_peak_is_gc = n_peak_is_gc
    )
}, mc.cores = env$N_CPU)
# })
names(sampleData) <- samples$sample_name

message("refactoring sample data types into matrices (row=bin x col=sample)")
n_obs_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, n_obs]
})
colnames(n_obs_bin_smp) <- samples$sample_name
mpp_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, mpp]
})
colnames(mpp_bin_smp) <- samples$sample_name
gc_bin_smp <- sapply(samples$sample_name, function(sample_name) {
    sampleData[[sample_name]]$bins[, gc]
})
colnames(gc_bin_smp) <- samples$sample_name
f_obs_ref_isl_smp <- list(primary = {
    m <- sapply(samples$sample_name, function(sample_name) {
        sampleData[[sample_name]]$f_obs_isl
    })
    colnames(m) <- samples$sample_name
    m
})
n_peak_is_gc_smp <- sapply(0:1, function(peak) {
    m1 <- sampleData[[1]]$n_peak_is_gc[[peak + 1]]
    array(
        do.call(c, lapply(samples$sample_name, function(sample_name) {
            as.vector(sampleData[[sample_name]]$n_peak_is_gc[[peak + 1]])
        })),
        dim = c(nrow(m1), ncol(m1), nSamples),
        dimnames = list(
            insert_size = NULL,
            gc_bin      = NULL,
            sample      = samples$sample_name
        )
    )
}, simplify = FALSE, USE.NAMES = TRUE)

message()
message("saving recollated sample data for app")
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
    samples = samples,
    references = readRDS(paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.'))$references,
    bins = bins[, .(
        chrom = chr, start0 = bs0, 
        included = included # where peak bins are now excluded (not included)
    )],
    # used by score.R
    n_obs_bin_smp = n_obs_bin_smp,
    n_wgt_bin_smp = NULL,
    mpp_bin_smp   = mpp_bin_smp,
    gc_bin_smp    = gc_bin_smp,
    f_obs_ref_isl_smp = f_obs_ref_isl_smp, # primary only
    # not used by score.R
    n_peak_is_gc_smp = n_peak_is_gc_smp
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "recollate.rds", sep = '.')
)
str(obj)
