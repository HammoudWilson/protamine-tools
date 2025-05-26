# action:
#     find positioned nucleosomes independently of TSSs (but score them based on proximity to TSSs)
#     logic is that positioned nucleosomes:
#         have high coverage by nucleosome-size fragments > 150bp (nuc)
#         have low coverage by subnucleosome-size fragments < 150bp (subnuc)
#         have few fragment endpoints of any length within them (ends)
#         therefore have a high nucleosome bias calculated as (nuc - subnuc) / (nuc + subnuc)
#         are highly accessible to Tn5 by virtue of flanking nucleosome-depleted regions
#         therefore have a high yield of nucleosome-size fragments relative to the rest of the genome
#      thus, positioned nucleosomes are identified as bins with:
#         high nucleosome bias
#         high nucleosome-length fragment coverage
#      and nucleosome-depleted regions are identified as bins with:
#         low (negative) nucleosome bias
#         high subnucleosome-length fragment coverage
# input:
#     output object from `find_nuc.R`
# outputs:
#     output object with positioned nucleosome information added as objects:
#         bins (list of scored bins for plotting in track)
#         calls (list of positioned nucleosomes and nucleosome-depleted regions)

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
        'ACTION_DIR',
        'DATA_FILE_PREFIX'
    ),
    integer = c(
        'MIN_MAPQ',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
BIN_SIZE <- 25L
MIN_NUC_SIZE <- 150L
MAX_NUC_SIZE <- 300L
BINS_PER_WINDOW  <- as.integer(MIN_NUC_SIZE / BIN_SIZE - 1L) # thus, 5-bin, 125bp windows
MIN_BASES_CALL    <- 10L * BIN_SIZE # base equivalent of 10 reads per bin
MIN_CALL_QUANTILE <- 0.9
MIN_CALL_BIAS     <- 0.9
#=====================================================================================

#=====================================================================================
# find positioned nucleosomes ab initio
#-------------------------------------------------------------------------------------

message("loading tssFrags object")
tssFragsFile <- paste(env$DATA_FILE_PREFIX, "tssFrags.rds", sep = '.')
tfd <- readRDS(tssFragsFile)
stages <- c('early_round', 'late_round') # unique(tfd$samples$stage)

message("scoring bins for nucleosome and subnucleosome base counts")
scoreAbInitio <- paste('bash', file.path(env$ACTION_DIR, 'score_ab_initio.sh'))
abInitioBins <- sapply(stages, function(stage_) {
    message(paste(' ', stage_))
    sample_names <- tfd$samples[stage == stage_, sample_name]
    bamFiles <- sapply(sample_names, function(sample_name_) {
        sample <- tfd$samples[sample_name == sample_name_]
        file.path(tfd$reference$input_dir, paste0(sample$filename_prefix, '.*.bam'))
    })
    chroms <- tfd$reference$chroms
    # chroms <- c("chr18","chr19")
    do.call(rbind, mclapply(chroms, function(chrom) { 
    # do.call(rbind, lapply(chroms, function(chrom) { 
        message(paste('   ', chrom))
        fread(cmd = paste(
            scoreAbInitio, 
            paste0('"', bamFiles, '"'), chrom, 
            BIN_SIZE, MIN_NUC_SIZE, MAX_NUC_SIZE, 
            tfd$reference$fai_file
        )) 
    }, mc.cores = env$N_CPU))
    # }))
}, simplify = FALSE, USE.NAMES = TRUE) 

message("finding positioned nucleosomes and nucleosome-depleted regions ab initio")
abInitioCalls <- sapply(stages, function(stage_) {
    message(paste(' ', stage_))
    chroms <- abInitioBins[[stage_]][, unique(chrom)]
    do.call(rbind, mclapply(chroms, function(chrom_) { 
    # do.call(rbind, lapply(chroms, function(chrom_) { 
        message(paste('   ', chrom_))
        bins <- abInitioBins[[stage_]][chrom == chrom_]
        bins[, ":="(
            isNucCall =      bias         >=  MIN_CALL_BIAS & 
                             nuc_mean     >  quantile(nuc_mean, MIN_CALL_QUANTILE) & 
                             nuc_mean     >= MIN_BASES_CALL,
            isSubNucCall =   bias         <= -MIN_CALL_BIAS & 
                             subnuc_mean  >  quantile(subnuc_mean, MIN_CALL_QUANTILE) & 
                             subnuc_mean  >= MIN_BASES_CALL
        )]
        if(any(bins$isNucCall) || any(bins$isSubNucCall)) {
            bins[isNucCall == TRUE | isSubNucCall == TRUE, # these are mutually exclusive, a bin can't be both (but can be neither)
                .(
                    chrom       = chrom_, 
                    start0      = min(start0), 
                    end1        = max(start0) + BIN_SIZE,
                    nuc_mean    = round(mean(nuc_bin),    1), 
                    subnuc_mean = round(mean(subnuc_bin), 1), 
                    bias        = round((sum(nuc_bin) - sum(subnuc_bin)) / (sum(nuc_bin) + sum(subnuc_bin)), 3)
                ), 
                keyby = rleid(isNucCall, isSubNucCall)
            ][, .(
                chrom, start0, end1, nuc_mean, subnuc_mean, bias
            )]
        } else {
            NULL
        }
    }, mc.cores = env$N_CPU))
    # }))
}, simplify = FALSE, USE.NAMES = TRUE) 

# TODO: use all fragments that flank each positioned nucleosome to finds its midpoint
# for now, can use the midpoint of the bins as a surrogate

# TODO: identify closely spaced calls and parse into concatenated call regions

message()
message("saving data for app")
ai <- list(
    binSize = BIN_SIZE,
    binsPerWindow = BINS_PER_WINDOW,
    windowCenterOffset = BINS_PER_WINDOW * BIN_SIZE / 2,
    binCenterOffset = BIN_SIZE / 2,
    bins = {
        x <- lapply(abInitioBins, function(x) x[, .(chrom, start0, nuc_bin, subnuc_bin)])
        names(x) <- names(abInitioBins)
        x
    },
    calls = abInitioCalls
)
saveRDS(
    ai, 
    file = paste(env$DATA_FILE_PREFIX, "abInitio.rds", sep = '.')
)
str(ai)
#=====================================================================================
