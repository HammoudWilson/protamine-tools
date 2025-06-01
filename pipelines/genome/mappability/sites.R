# action:
#   count all possible discontiguous Tn5 nonamers within genome mappable regions at different insert size levels
# input:
#   mappability files created by mappability.sh
# outputs:
#   for use by atac/collate when assembling bin counts:
#       ${GENOME_METADATA_PREFIX}.<refType>.tn5_site_freqs_exp.txt.gz
#           one file per reference genome (refType = primary, spike_in)
#           table of expected Tn5 nonamer frequencies by insert size level
#   for use in app:
#       ${GENOME_METADATA_PREFIX}.<refType>.tn5_site_counts_exp.rds
#           one file per reference genome (refType = primary, spike_in)
#           complete view of Tn5 nonamer counts by insert size level
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
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'GENOME_FASTA_SHM',
        'MAPPABILITY_SIZE_LEVELS',
        'SIZE_LEVEL_MAXIMA',
        'MAPPABILITY_FILE_PREFIX',
        'TN5_PREFERENCE_POSITIONS',
        'TN5_KMERS_FILE',
        'ACTION_DIR',
        'TASK_DIR',
        'GENOME_METADATA_PREFIX'
    ),
    integer = c(
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
# loop through all insert size level to extract insert cleavage sites
#-------------------------------------------------------------------------------------
insertSizeLevels <- as.integer(strsplit(env$MAPPABILITY_SIZE_LEVELS, " ")[[1]])
sizeLevelMaxima  <- as.integer(strsplit(env$SIZE_LEVEL_MAXIMA, " ")[[1]])
nSizeLevels <- length(insertSizeLevels)

message("counting possible genome cleavage sites by insert size level")
collectGenomeSites <- paste('bash', file.path(env$ACTION_DIR, 'collect_genome_sites.sh'))
n_tn5 <- mclapply(1:nSizeLevels, function(sizeLevelI) {
# n_tn5 <- lapply(1:nSizeLevels, function(sizeLevelI) {
    insertSizeLevel <- insertSizeLevels[sizeLevelI]
    maxInsertSize <- sizeLevelMaxima[sizeLevelI]
    label <- paste0('   insert size level = ', insertSizeLevel)
    message(label)
    mppFile <- paste0(env$MAPPABILITY_FILE_PREFIX, ".k_", insertSizeLevel, ".e_*.bed.gz")
    n_tn5 <- fread(
        cmd = paste(collectGenomeSites, mppFile, insertSizeLevel, maxInsertSize), 
        sep = "\t", header = TRUE
    )
    # this loop only calculates one insert size level at a time, but collect_cleavage_sites.sh returns all size levels...
# Classes ‘data.table’ and 'data.frame':  524288 obs. of  21 variables:
#  $ genome : chr  "mm39" "mm39" "mm39" "mm39" ...
#  $ tn5_35 : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ tn5_40 : int  1 0 0 0 0 0 0 0 0 0 ...
    tn5_N <- paste0("tn5_", insertSizeLevel)
    list(
        # each genome's has yields a vector of nonamer counts for this insert size level
        primary  = n_tn5[genome == env$PRIMARY_GENOME,  ..tn5_N][[1]],
        spike_in = n_tn5[genome == env$SPIKE_IN_GENOME, ..tn5_N][[1]]
    )
}, mc.cores = env$N_CPU)
# })
#=====================================================================================

#=====================================================================================
# parse to final outputs
#-------------------------------------------------------------------------------------
nSites <- length(n_tn5[[1]]$primary)
for(refType in c("primary", "spike_in")) {
    message()
    filename <- "tn5_site_counts_exp.rds"
    message(paste("refactoring ", refType, "to", filename))
    n_tn5_m <- sapply(1:nSizeLevels, function(sizeLevelI) {
        n_tn5[[sizeLevelI]][[refType]]
    })
    saveRDS(
        n_tn5_m, 
        file = paste(env$GENOME_METADATA_PREFIX, refType, filename, sep = '.')
    )
    str(n_tn5_m)

    filename <- "tn5_site_freqs_exp.txt.gz"
    message(paste("refactoring ", refType, "to", filename))
    f_tn5_m <- sapply(1:nSizeLevels, function(sizeLevelI) {
        n_tn5_sum <- sum(n_tn5_m[, sizeLevelI])
        if(n_tn5_sum > 0) {
            n_tn5_m[, sizeLevelI] / n_tn5_sum
        } else {
            rep(0, nSites)
        }
    })
    fwrite(
        as.data.table(f_tn5_m), 
        file = paste(env$GENOME_METADATA_PREFIX, refType, filename, sep = '.'), 
        # forgot to set sep = "\t" in fwrite, so it defaults to comma, i.e., a CSV file
        quote = FALSE,
        row.names = FALSE, 
        col.names = FALSE
    )
    str(f_tn5_m)
}
#=====================================================================================
