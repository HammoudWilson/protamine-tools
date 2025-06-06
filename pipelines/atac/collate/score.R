# action:
#     calculate and analyze the distributions of the following score metrics per bin:
#         genome-level (not specific to samples at different spermiogenic stages):
#               gc       = aggregated bin GC content over all samples
#               txn      = (nascent) Transcription count per million reads
#         sample-level (specific to each sample at different spermiogenic stages):
#               gcrz_obs = GC Residual Z-Score (excess or deficit of reads relative to GC peers), based on observed read counts
#               gcrz_wgt = as above, now based on Tn5 site weights, i.e., including Tn5 kmer preference correction before GC regression
#               nrll     = Normalized Relative Log Likelihood of protamine- vs. histone-associated insert sizes
#     scores are only calculated relative to the primary genome distribution
#         bin-level outputs are adjusted to only include primary genome bins
#         the bins chrom field is stripped of the genome name for indexing by browser
# input:
#     sample metadata file
#     composite genome bins BED file
#     insert_spans file created by atac/sites
#     output of the atac/collate/collate action step
# outputs:
#     a gcBiasModels RDS file with GC bias regression models for each sample 
#     a scores RDS file including a nested list with metadata about all scores aggregated by:
#         individual sample
#         spermiogenic stage
#         the difference between early and late spermiogenic stages
#         where each level has:
#             the distribution of bin scores in all included autosomal bins
#             the peak, median, mean, and standard deviation of the distribution
#     a series of tabixed files for each score type, which can be used in the app browser
#         one file per score type, per data type
#         with columns for each sample, aggregated stage, and aggregated stageType
#         score types are as above
#         data types are raw score values, z-scores, and quantiles across included autosomal bins
#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
    library(MASS) # for glm.nb
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
        'DATA_NAME',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE',
        'STAGE_TYPES',
        'GC_LIMITS',
        'TRANSCRIPTION_BED',
        'SHM_FILE_PREFIX',
        'TMP_FILE_PREFIX'
    ),
    integer = c(
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'BIN_SIZE',
        'N_CPU'
    )
))
if(env$TRANSCRIPTION_BED == "NA") env$TRANSCRIPTION_BED <- paste(env$DATA_FILE_PREFIX, "nascent_transcriptome_unstranded.bed.gz", sep = '.')
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$MODULES_DIR, 'score')
sourceScripts(rUtilDir, c('score_functions'))
rUtilDir <- file.path(env$MODULES_DIR, 'score', 'nbinomCountsGC')
sourceScripts(rUtilDir, c('nbinomCountsGC_class', 'nbinomCountsGC_methods'))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# loop through all BAM files to determine their coverage in each genome bin
#-------------------------------------------------------------------------------------
message("loading collate step output")
collate <- readRDS(paste(env$DATA_FILE_PREFIX, "collate.rds", sep = '.'))

message("parsing spermiogenic stage types and GC limits")
stageTypes <- unpackStageTypes(env)
gcLimits <- strsplit(env$GC_LIMITS, ",")[[1]]

message("restricting attention to primary genome bins")
isPrimaryGenome <- collate$bins[, endsWith(chrom, paste0("-", env$PRIMARY_GENOME))]

message("calculating aggregated bin GC content and number of alleles")
collate$bins[isPrimaryGenome, gc := rowMeans(collate$gc_bin_smp[isPrimaryGenome,], na.rm = TRUE)]
collate$bins[isPrimaryGenome, nAlleles := ifelse(chrom %in% paste0(c('chrX', 'chrY'), "-", env$PRIMARY_GENOME), 1L, 2L)]

message("restricting attention to included primary autosomal bins for score calculations")
isIncluded     <- collate$bins[, included == 1]
isAutosome     <- collate$bins[, !startsWith(chrom, c("chrX-", "chrY-"))]
passesGcLimits <- collate$bins[, gc >= gcLimits[1] & gc <= gcLimits[2]]
includedAutosomeBins <- isPrimaryGenome & isIncluded & isAutosome & passesGcLimits

message("extracting histone- and protamine-associated insert size distributions for primary genome")
emissProbsFile <- extractInsertSizeEps(collate$samples, collate$f_obs_ref_isl_smp$primary, env)

message("analyzing genome-level scores")
# add columns to collate$bins for genome-level scores for use in browser (written below to data package)
scores <- list(genome = list())
message("  gc")
scores$genome$gc <- analyzeScoreDist(collate$bins$gc, scoreTypes$genome$gc, scoreTypes$genome$gc$name)
collate$bins[isPrimaryGenome, gc_z := readRDS(scores$genome$gc$scoreFile)$z]
scores$genome$txn <- if(file.exists(env$TRANSCRIPTION_BED)) {
    message(paste("  txn:", env$TRANSCRIPTION_BED))
    txn_cpm <- fread(env$TRANSCRIPTION_BED)[[4]] # requires BED4 with normalized bin transcription value in column 4
    txn <- analyzeScoreDist(txn_cpm, scoreTypes$genome$txn, scoreTypes$genome$txn$name, return_scores = TRUE)
    collate$bins[isPrimaryGenome, txn := txn$scores[isPrimaryGenome]] # txn in bins is log10_cpm
    txn$scores <- NULL
    txn
} else {
    message("  no txn data")
    collate$bins[isPrimaryGenome, txn := NA_real_]
    NULL
}

message("analyzing sample-level scores")
scores$sample <- sapply(names(scoreTypes$sample), function(scoreTypeName) {
    scoreType <- scoreTypes$sample[[scoreTypeName]]
    # if(scoreType$gcBiasDependent) return(NULL)
    message(paste(" ", scoreTypeName))
    scoreFnName <- paste("get", scoreTypeName, sep = "_")
    sampleScores <- analyzeSampleScores(get(scoreFnName), scoreType, emissProbsFile)
    aggregateScores <- aggregateSampleScores(sampleScores, scoreType)
    # now that aggregation is complete, remove scores from sampleScores
    for(sample_name in names(sampleScores)) sampleScores[[sample_name]]$scores <- NULL
    list(
        sampleScores    = sampleScores,
        aggregateScores = aggregateScores
    )
}, simplify = FALSE, USE.NAMES = TRUE)

message("aggregating and saving GC regression fits for app")
gcBiasModels <- sapply(c("gcrz_obs", "gcrz_wgt"), function(grcz_type) {
    sapply(collate$samples$sample_name, function(sample_name) {
        tmpFile <- paste(env$SHM_FILE_PREFIX, "gcBiasModel", grcz_type, sample_name, "rds", sep = '.')
        readRDS(tmpFile)
    }, simplify = FALSE, USE.NAMES = TRUE)
}, simplify = FALSE, USE.NAMES = TRUE)
saveRDS(
    gcBiasModels, 
    file = paste(env$DATA_FILE_PREFIX, "gcBiasModels.rds", sep = '.')
)

message("refactoring bin-level score data for tabixed retrieval in app")
env$SCORES_DIR <- file.path(env$TASK_DIR, "scores")
scoresPrefix <- file.path(env$SCORES_DIR, env$DATA_NAME)
if(!dir.exists(env$SCORES_DIR)) dir.create(env$SCORES_DIR)
templateDt <- collate$bins[isPrimaryGenome, .(
    chrom  = sub(paste0("-", env$PRIMARY_GENOME), "", chrom),
    start0 = start0, 
    end1   = start0 + env$BIN_SIZE
)]
scoreTables <- list()
for(scoreTypeName in names(scoreTypes$sample)){
    scoreType <- scoreTypes$sample[[scoreTypeName]]
    # if(scoreType$gcBiasDependent) return(NULL)
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
            fill_dt(scores$sample[[scoreTypeName]]$aggregateScores$by_stageType)
            fill_dt(scores$sample[[scoreTypeName]]$aggregateScores$delta)
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
message("saving output for app") # see above for saveRDS of gcBiasModels
env$MAPPABILITY_SIZE_LEVELS <- as.integer(strsplit(env$MAPPABILITY_SIZE_LEVELS, " ")[[1]])
obj <- list(
    env = env[c(
        'PRIMARY_GENOME',
        'SPIKE_IN_GENOME',
        'GENOME',
        'MAPPABILITY_SIZE_LEVELS',
        'MIN_INSERT_SIZE',
        'MAX_INSERT_SIZE',
        'BIN_SIZE',
        'HISTONE_STAGE',
        'PROTAMINE_STAGE',
        'SCORES_DIR'
    )],
    samples     = collate$samples,
    references  = collate$references,
    stageTypes  = stageTypes,
    gcLimits    = gcLimits,
    scores      = scores, # summary statistics only, see bed.bgz files for sample scores, bins for genome scores
    scoreTables = scoreTables,
    bins = collate$bins[isPrimaryGenome, .(
        chrom    = sub(paste0("-", env$PRIMARY_GENOME), "", chrom),
        start0   = start0,
        included = included,
        gc       = gc,
        gc_z     = gc_z,
        txn      = txn
    )]
)
saveRDS(
    obj, 
    file = paste(env$DATA_FILE_PREFIX, "scores.rds", sep = '.')
)
str(obj)
#=====================================================================================
