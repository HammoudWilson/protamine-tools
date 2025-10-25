# action:
#     explore the spatio-temporal relationship of positioned nucleosomes found ab initio
#     the goal is to establish localization (anti)patterns of early and late accessibility peaks
# input:
#     ab_initio.rds from ab_initio.R
# outputs:
#     lag_analysis.rds for use in the app

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
        'DATA_FILE_PREFIX',
        'TASK_DIR',
        'MDI_DIR'
    ),
    integer = c(
        'N_CPU'
    )
))
env$MAX_LAG <- 10e6
env$MIN_LAG_PAIRS <- 100
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# lag analysis
#-------------------------------------------------------------------------------------
message("loading ab initio data")
aiFile <- paste(env$DATA_FILE_PREFIX, "ab_initio.rds", sep = '.')
ai <- readRDS(aiFile)

message("parsing stages, early_RS to early_ES only")
stage_rpkm_cols   <- paste(ai$stages[1:5], "rpkm",   sep = '_')
stage_scaled_cols <- paste(ai$stages[1:5], "scaled", sep = '_')

message("extracting required data = overlap groups that passed RPKM filters")
regions <- ai$regions[index_stage == 'overlap_group' & passed_rpkm_filters == TRUE]
rm(ai)

message("scaling regions as log2 fold-change from the mean of all stages") # row-wise normalization
scaled <- do.call(rbind, mclapply(1:nrow(regions), function(i) {
    rpkm <- unlist(regions[i, ..stage_rpkm_cols])
    mean_rpkm <- mean(rpkm, na.rm = TRUE)
    data.table(t(log2(rpkm / mean_rpkm))) # can yield -Inf if stage rpkm == 0 in a region; must be handled downstream; mean_rpkm cannot be 0
}, mc.cores = env$N_CPU))
setnames(scaled, stage_scaled_cols)
regions <- cbind(regions[, .(chrom, start0, end1, stage_mean)], scaled)
rm(scaled)

message("calculating peak lags")
regions[, center := (start0 + end1) / 2]
lag_data <- do.call(rbind, lapply(unique(regions$chrom), function(chrom_){
    message(paste("  ", chrom_))
    cint <- regions[chrom == chrom_][order(center)]
    if(nrow(cint) < 2) return(NULL)
    stageData <- as.matrix(cint[, .SD, .SDcols = stage_scaled_cols])
    do.call(rbind, mclapply(1:(nrow(cint) - 1), function(i){
        j <- which(cint$center > cint$center[i] & cint$center < cint$center[i] + env$MAX_LAG)
        if(length(j) == 0) return(NULL)
        data.table(
            log10dist = as.integer(round(log10(cint$center[j] - cint$center[i]) * 10, 0)), # so, e.g., 1kb => 30, etc.
            diff2 = (cint$stage_mean[j] - cint$stage_mean[i])^2,
            corr = sapply(j, function(jj) cor(stageData[i, ], stageData[jj, ])),
            diff = abs(cint$stage_mean[j] - cint$stage_mean[i])
        )
    }, mc.cores = env$N_CPU))
}))[, .(
        mean_corr = mean(corr, na.rm = TRUE),
        sd_corr = sd(corr, na.rm = TRUE),
        var = sum(diff2) / .N / 2,
        mean_diff = mean(diff, na.rm = TRUE),
        sd_diff = sd(diff, na.rm = TRUE),
        n_pairs = .N
    ), 
    keyby = .(log10dist)
][
    n_pairs >= env$MIN_LAG_PAIRS
]

message()
message("saving data for app")
lag_data[, log10dist := log10dist / 10] # back to log10(bp)
saveRDS(
    lag_data, 
    file = paste(env$DATA_FILE_PREFIX, "lag_analysis.rds", sep = '.')
)
str(lag_data)
#=====================================================================================
