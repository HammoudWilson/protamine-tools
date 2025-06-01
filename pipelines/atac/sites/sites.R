# action:
#   count the discontiguous Tn5 nonamers at all observed insert cleavage sites in parallel over multiple samples
# input:
#   sample metadata file
#   aligned, sorted, and indexed bam files
#       one per sample
#       aligned to composite genome
#       not deduplicated, we will dedup here
# outputs:
#   for use by atac/collate when assembling bin counts:
#       ${TASK_DIR}/insert_spans/*.bed.gz
#           one file per sample
#           restricted to fully filtered insert spans
#           BED3 + insert_size_level + left_tn5_nonamer + right_tn5_nonamer
#       ${DATA_FILE_PREFIX}.<refType>.tn5_site_freqs_obs.txt.gz
#           one file per reference genome (refType = primary, spike_in)
#           table of observed Tn5 nonamer frequencies aggregated over all samples by insert size level
#   for use in app:
#       ${DATA_FILE_PREFIX}.<refType>.tn5_site_counts_obs.rds
#           one file per reference genome (refType = primary, spike_in)
#           complete view of Tn5 nonamer counts by insert size level and sample
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
        'GENOME_GAPS_FILE',
        'GENOME_EXCLUSIONS_BED',
        'INSERT_SPANS_DIR',
        'MAPPABILITY_SIZE_LEVELS',
        'MAPPABILITY_FILE_PREFIX',
        'TN5_PREFERENCE_POSITIONS',
        'TN5_KMERS_FILE',
        'ACTION_DIR',
        'TMP_FILE_PREFIX',
        'DATA_FILE_PREFIX',
        'TASK_DIR' 
    ),
    integer = c(
        'MIN_MAPQ',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'TN5_FLANK_PADDING',
        'N_TN5_PREFERENCE_POSITIONS',
        'TN5_PREFERENCE_SPAN',
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to extract insert cleavage sites
#-------------------------------------------------------------------------------------
message("loading and sorting sample metadata")
samples <- fread(env$METADATA_FILE)[order(staging_order)]
nSamples <- nrow(samples)

message("collecting insert cleavage sites by sample")
collectCleavageSites <- paste('bash', file.path(env$ACTION_DIR, 'collect_cleavage_sites.sh'))
n_tn5 <- mclapply(1:nrow(samples), function(sampleI) {
# n_tn5 <- lapply(1:nrow(samples), function(sampleI) {
    smp <- samples[sampleI]
    label <- paste0('   ', smp$filename_prefix, ' = ', smp$sample_name, " (", smp$staging, ")")
    message(label)
    bamFile <- file.path(env$INPUT_DIR, paste0(smp$filename_prefix, '.*.bam'))
    tmpPrefix <- paste0(env$TMP_FILE_PREFIX, ".atac_sites.", smp$filename_prefix)
    n_tn5 <- fread(
        cmd = paste(collectCleavageSites, bamFile, tmpPrefix, smp$sample_name, smp$filename_prefix), 
        sep = "\t", header = TRUE
    )
# Classes ‘data.table’ and 'data.frame':  524288 obs. of  21 variables:
#  $ genome : chr  "mm39" "mm39" "mm39" "mm39" ...
#  $ tn5_35 : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ tn5_40 : int  1 0 0 0 0 0 0 0 0 0 ...
    # collect insert counts by stage
    n_ins_stage <- fread(paste0(tmpPrefix, ".n_ins_stage.txt"))
    setnames(n_ins_stage, c('sample_name', 'stage', 'n'))
    list(
        # each genome's matrix has rows = nonamer, cols = insert_size_level
        primary  = as.matrix(n_tn5[genome == env$PRIMARY_GENOME,  -1]),
        spike_in = as.matrix(n_tn5[genome == env$SPIKE_IN_GENOME, -1]),
        n_ins_stage = n_ins_stage
    )
}, mc.cores = env$N_CPU)
# })
names(n_tn5) <- samples$sample_name
#=====================================================================================

#=====================================================================================
# parse to final outputs
#-------------------------------------------------------------------------------------
message()
message("writing summary counts by filtering stage")
summary <- do.call(rbind, lapply(samples$sample_name, function(sample_name) {
    n_tn5[[sample_name]]$n_ins_stage
}))
summary <- dcast(summary, sample_name ~ stage, value.var = 'n', fun.aggregate = sum, fill = 0L)
summary[, ":="(
    frac_included_2_1 = round(stage_2_included / stage_1_proper_pairs, 3),
    frac_unique_3_2   = round(stage_3_dedup    / stage_2_included, 3),
    frac_mappable_4_3 = round(stage_4_mappable / stage_3_dedup, 3)
)]
fwrite(
    summary, 
    file = paste(env$DATA_FILE_PREFIX, "n_inserts_by_stage.txt", sep = '.'), 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
)

nSites      <- nrow(n_tn5[[1]]$primary)
nSizeLevels <- ncol(n_tn5[[1]]$primary)
for(refType in c("primary", "spike_in")) {
    message()
    filename <- "tn5_site_counts_obs.rds"
    message(paste("refactoring ", refType, "to", filename))
    n_tn5_arr <- array(
        do.call(c, lapply(samples$sample_name, function(sample_name) {
            as.vector(n_tn5[[sample_name]][[refType]])
        })),
        dim = c(
            nSites, 
            nSizeLevels, 
            nSamples
        ),
        dimnames = list(
            tn5_site = NULL,
            insert_size_level = NULL,
            sample = samples$sample_name
        )
    )
    saveRDS(
        n_tn5_arr, 
        file = paste(env$DATA_FILE_PREFIX, refType, filename, sep = '.')
    )
    str(n_tn5_arr)

    filename <- "tn5_site_freqs_obs.txt.gz"
    message(paste("refactoring ", refType, "to", filename))
    f_tn5_all <- sapply(1:nSizeLevels, function(sizeLevelI) {
        n_tn5_all <- rowSums(n_tn5_arr[, sizeLevelI, ])
        n_tn5_sum <- sum(n_tn5_all)
        if(n_tn5_sum > 0) {
            n_tn5_all / n_tn5_sum
        } else {
            rep(0, nSites)
        }
    })
    fwrite(
        as.data.table(f_tn5_all), 
        file = paste(env$DATA_FILE_PREFIX, refType, filename, sep = '.'), 
        # forgot to set sep = "\t" in fwrite, so it defaults to comma, i.e., a CSV file
        quote = FALSE,
        row.names = FALSE, 
        col.names = FALSE
    )
    str(f_tn5_all)
}
#=====================================================================================
