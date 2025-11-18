# action:
#     calculate and analyze the distributions of the following score metrics per bin:
#         sample-level (specific to each sample at different spermiogenic stages):
#               CUt&Tag RPKM for different antibody targets
#     scores are only calculated relative to the primary genome distribution
#         bin-level outputs are adjusted to only include primary genome bins
#         the bins chrom field is stripped of the genome name for indexing by browser
# input:
#     Cut&Tag sample metadata file
#     output of `atac collate`, in particular the genome bins file
#     Cut&Tag bigwig count files compiled outside of this pipeline (see example command in score_functions.R)
# outputs:
#     a scores RDS file including a nested list with metadata about all scores aggregated by:
#         individual sample
#         spermiogenic stage
#         where each level has:
#             bin scores as unnormalized read counts (normalization is done in the app)
#             the distribution of bin scores in all included autosomal bins
#             the peak, median, mean, and standard deviation of the distribution
#     a series of tabixed files for each score type, which can be used in the app browser
#         one file per score type, per data type
#         with columns for each sample and aggregated stage
#         data types are raw score values and quantiles across included autosomal bins
#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(rtracklayer)
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
        'GENOME_BINS_BED',
        'ACTION_DIR',
        'DATA_FILE_PREFIX',
        'DATA_NAME',
        'GC_LIMITS',
        'SHM_FILE_PREFIX',
        'TMP_FILE_PREFIX'
    ),
    integer = c(
        'BIN_SIZE',
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$MODULES_DIR, 'score')
sourceScripts(rUtilDir, c('score_functions'))
sourceScripts(env$ACTION_DIR, c('score_functions')) # additional and overrides specific to Cut&Tag data
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------
message("loading collate step output")
collateFile <- paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.')
collate <- readRDS(collateFile)
bins <- collate$bins
gc_bin_smp <- collate$gc_bin_smp
rm(collate)

message("parsing spermiogenic stages")
gcLimits <- strsplit(env$GC_LIMITS, ",")[[1]]

message("restricting attention to primary genome bins")
isPrimaryGenome <- bins[, endsWith(chrom, paste0("-", env$PRIMARY_GENOME))]

message("calculating aggregated number of alleles")
bins[isPrimaryGenome, gc := rowMeans(gc_bin_smp[isPrimaryGenome,], na.rm = TRUE)]
bins[isPrimaryGenome, nAlleles := ifelse(chrom %in% paste0(c('chrX', 'chrY'), "-", env$PRIMARY_GENOME), 1L, 2L)]

message("restricting attention to included primary autosomal bins for score calculations")
isIncluded     <- bins[, included == 1]
isAutosome     <- bins[, !startsWith(chrom, c("chrX-", "chrY-"))]
passesGcLimits <- bins[, gc >= gcLimits[1] & gc <= gcLimits[2]]
includedAutosomeBins <- isPrimaryGenome & isIncluded & isAutosome & passesGcLimits
primaryIncludedBins <- isPrimaryGenome & isIncluded

message("loading Cut&Tag sample metadata")
cuttagSamples <- fread(env$METADATA_FILE, sep = ",", header = TRUE)

message("analyzing sample-level scores")
scores <- list() # maintain the same file structure as ATAC sample scores
scores$sample <- sapply(names(scoreTypes$sample), function(antibody_target){
    message(paste(" ", antibody_target))
    scoreType <- scoreTypes$sample[[antibody_target]]
    sampleScores <- analyzeSampleScores(scoreType, antibody_target)
    aggregateScores <- aggregateSampleScores(sampleScores, scoreType, antibody_target)
    # now that aggregation is complete, remove scores from sampleScores
    for(sample_name in names(sampleScores)) sampleScores[[sample_name]]$scores <- NULL
    list(
        sampleScores    = sampleScores,
        aggregateScores = aggregateScores
    )
}, simplify = FALSE, USE.NAMES = TRUE)

message("refactoring bin-level score data for tabixed retrieval in app")
env$SCORES_DIR <- file.path(env$TASK_DIR, "cuttag")
scoresPrefix <- file.path(env$SCORES_DIR, env$DATA_NAME)
if(!dir.exists(env$SCORES_DIR)) dir.create(env$SCORES_DIR)
templateDt <- bins[isPrimaryGenome, .(
    chrom  = sub(paste0("-", env$PRIMARY_GENOME), "", chrom),
    start0 = start0, 
    end1   = start0 + env$BIN_SIZE
)]
scoreTables <- list()
for(scoreTypeName in names(scoreTypes$sample)){
    scoreType <- scoreTypes$sample[[scoreTypeName]]
    if(!is.list(scores$sample[[scoreTypeName]])){
        next
    }
    message(paste(" ", scoreTypeName))
    scoreTables[[scoreTypeName]] <- list()
    for (dataType in c("score","z","quantile")) {
        if(dataType == "score" || dataType %in% scoreType$include){
            dt <- copy(templateDt)
            colNames <- names(dt)
            fill_dt <- function(wrkScores){
                for (dataSet in names(wrkScores)) dt[[dataSet]] <<- readRDS(wrkScores[[dataSet]]$scoreFile)[[dataType]]
                colNames <<- c(colNames, names(wrkScores))
            }
            fill_dt(scores$sample[[scoreTypeName]]$sampleScores)
            fill_dt(scores$sample[[scoreTypeName]]$aggregateScores$by_stage)
            gzFile  <- paste(env$SHM_FILE_PREFIX, "score_table_tmp.gz", sep = '.')
            bgzFile <- paste(scoresPrefix, scoreTypeName, dataType, "bed.bgz", sep = '.')
            fwrite(dt, file = gzFile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
            system(sprintf('zcat %s | bgzip -c > %s', gzFile, bgzFile))
            system2("tabix", args = c("-p", "bed", bgzFile))
            unlink(gzFile)
            scoreTables[[scoreTypeName]][[dataType]] <- list(
                bgzFile = basename(bgzFile),
                colNames = colNames
            )
        }
    }
}

message()
message("saving output for app")
obj <- list(
    # analogous to collate.rds
    env = env[c(
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'BIN_SIZE',
        'SCORES_DIR'
    )],
    samples     = cuttagSamples,
    # analogous to score.rds
    gcLimits    = gcLimits,
    scores      = scores, # summary statistics only, see bed.bgz files for sample scores, bins for genome scores
    scoreTables = scoreTables
)
saveRDS(
    obj, 
    file = paste0(env$DATA_FILE_PREFIX, ".cuttag.rds")
)
str(obj)
#=====================================================================================
