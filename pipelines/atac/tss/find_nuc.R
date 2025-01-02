# action:
#     attempt to find positioned nucleosomes distal to active TSSs
#     require sufficient mononucleosome fragments in a focused region to establish its median position
# input:
#     output object from `tss.R`
# outputs:
#     updated TSS output object with positioned nucleosome information added to tss and tssFrags
#     only some active TSSs will have positioned nucleosomes identified and located, others are NA

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
        'DATA_FILE_PREFIX'
    ),
    integer = c(
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
MIN_NUC_FRAGS <- 25
# limits of positioned mononucleosomes downstream of TSS were deterimined empirically
MIN_NUC_TSS_DIST <- 0   # x-axis on V plots
MAX_NUC_TSS_DIST <- 250 
MIN_NUC_SIZE <- 150     # y-axis on V plots
MAX_NUC_SIZE <- 275
# the DIAG line prevents fragments in large nucleosome-free free regions from being considered
# effecttively a simple way to limit to reads in a rectangular region with the lower left corner cut off
DIAG_X1 <- 0 
DIAG_Y1 <- 210
DIAG_X2 <- 135
DIAG_Y2 <- 125
#-------------------------------------------------------------------------------------
# check if points fall below a line defined by two points
is_above_line <- function(x, y, x1, x2, y1, y2) {
    slope <- (y2 - y1) / (x2 - x1)
    intercept <- y1 - slope * x1
    y > slope * x + intercept
}
#=====================================================================================

#=====================================================================================
# find positioned nucleosomes
#-------------------------------------------------------------------------------------

message("loading tssFrags object")
tssFragsFile <- paste(env$DATA_FILE_PREFIX, "tssFrags.rds", sep = '.')
tfd <- readRDS(tssFragsFile)

message("setting TSS indices") # didn't do this in tss.R, used here and in app
tfd$tss$active$tss_i1   <- 1:nrow(tfd$tss$active)
tfd$tss$inactive$tss_i1 <- 1:nrow(tfd$tss$inactive)

message("locating positioned nucleosomes distal to early round active TSSs")
earlyRoundSamples <- tfd$samples[stage == "early_round", sample_name]
earlyRoundFrags <- do.call(rbind, tfd$tssFrags$active[earlyRoundSamples])
tfd$tss$active$nucleosome_distance <- unlist(mclapply(tfd$tss$active$tss_i1, function(tss_i1_) {
    frags <- earlyRoundFrags[tss_i1 == tss_i1_]
    if(nrow(frags) < MIN_NUC_FRAGS) return(NA_integer_)
    midpoint_tss <- frags[, .(
        midpoint_tss = start_to_tss + (end_to_tss - start_to_tss) / 2, 
        size = end1 - start0
    )][
        midpoint_tss %between% c(MIN_NUC_TSS_DIST, MAX_NUC_TSS_DIST) &
        size         %between% c(MIN_NUC_SIZE,     MAX_NUC_SIZE) &
        is_above_line(midpoint_tss, size, DIAG_X1, DIAG_X2, DIAG_Y1, DIAG_Y2),
        midpoint_tss
    ]
    if(length(midpoint_tss) < MIN_NUC_FRAGS) NA_integer_ 
    else as.integer(median(midpoint_tss))
}, mc.cores = env$N_CPU))

nActiveTss <- nrow(tfd$tss$active)
nPosNuc <- sum(!is.na(tfd$tss$active$nucleosome_distance))
meanNucDist <- mean(tfd$tss$active$nucleosome_distance, na.rm = TRUE)
sdNucDist   <-   sd(tfd$tss$active$nucleosome_distance, na.rm = TRUE)
message(paste(" found positioned nucleosomes for", nPosNuc, "of", nActiveTss, "active TSSs"))
message(paste(" mean +/- sd2 distance to positioned nucleosome:", round(meanNucDist, 0), "+/-", round(sdNucDist * 2, 0)))

message()
message("saving TSS inserts for app")
saveRDS(
    tfd, 
    file = tssFragsFile
)
str(tfd)
#=====================================================================================
