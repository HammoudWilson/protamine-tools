# action:
#   xxxxx
# input:
#     aligned and sorted bam files
# outputs:
#     xxxxx

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
        'INPUT_DIR',
        'GENOME_BINS_BED',
        'GENOME'
    ),
    integer = c(
        'BIN_SIZE'
    )
))
# #-------------------------------------------------------------------------------------
# # source R scripts
# rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
# sourceScripts(rUtilDir, c('utilities'))
# sourceScripts(file.path(rUtilDir, 'genome'), c('chroms', 'stats'))
# sourceScripts(file.path(rUtilDir, 'sequence'), c('general'))
# setCanonicalChroms()
# nonNuclear   <- c("chrM", "chrEBV")
# nonAutosomes <- c("chrX", "chrY", nonNuclear)
# nuclearChroms <- canonicalChroms[!(canonicalChroms %in% nonNuclear)]
# autosomes <- canonicalChroms[!(canonicalChroms %in% nonAutosomes)]
# env$MAX_INSERT_SIZE <- getMaxInsertSize()
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# feedback functions
#-------------------------------------------------------------------------------------

#=====================================================================================