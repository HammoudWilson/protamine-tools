# utilities for parsing bin scores and distributions

# sample-level score type functions
replaceNaN <- function(x){
    x[is.nan(x)] <- NA
    x
}
get_cpm <- function(bd, sample_name, ...){
    binCounts <- rowSums(bd$binCounts$genome[,, sample_name], na.rm = TRUE)
    replaceNaN(binCounts / sum(binCounts, na.rm = TRUE) * 1e6)
}
get_snif <- function(bd, sample_name, ...){
    replaceNaN(
        bd$binCounts$genome[,"subnucleosomal", sample_name] / 
        rowSums(bd$binCounts$genome[,, sample_name], na.rm = TRUE)
    )
}
get_nrll <- function(bd, sample_name_, emissProbsFile){
    ref <- bd$references$genome
    filename_prefix <- bd$samples[sample_name == sample_name_, filename_prefix]
    bamFile <- file.path(ref$input_dir, paste0(filename_prefix, '.*.bam'))
    script <- paste('bash', file.path(env$ACTION_DIR, 'get_bin_NRLL.sh'))
    unlist(lapply(ref$chroms, function(chrom) { # all bins values over all ordered chroms
        fread(cmd = paste(script, bamFile, chrom, ref$fai_file, emissProbsFile)) # one value per bin on chrom
    }))
}

# score types metadata
scoreTypes <- list(
    genome = list(
        gc = list(
            name = "Bin Fraction GC",
            gbBiasDependent = FALSE,
            distUnit = 0.01,
            class = "baseComposition"
        ),
        rpkm = list(
            name = "Transcription RPKM",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "transcription"
        )
    ),
    sample = list(
        cpm = list(
            name = "Counts Per Million",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "coverage"
        ),
        gcrz = list(
            name = "GC-Residual Z-Score",
            gbBiasDependent = TRUE, # thus, cannot be assessed until GC bias is established in app
            distUnit = 0.1,
            class = "coverage"
        ),
        snif = list(
            name = "Subnucleosomal Insert Fraction",
            gbBiasDependent = FALSE,
            distUnit = 0.01,
            class = "insertSize"
        ),
        nrll = list(
            name = "Subnucleosomal vs. Nucleosomal NRLL",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "insertSize"
        )
    )
)

# unpack the stage types from the user option value
unpackStageTypes <- function(env){
    stageTypes <- list()
    for(stageType in strsplit(env$STAGE_TYPES, ';')[[1]]){
        x <- strsplit(stageType, ':')[[1]]
        stageTypes[[x[1]]] <- strsplit(x[2], ',')[[1]]
    }
    stageTypes
}

# extract the nucleosomal and subnucleosomal insert size distributions
getStateEmissProbs <- function(isd, stage_){ # where a single specific stage is taken as being a sufficiently pure representation of a state
    stage_samples <- isd$samples[stage == stage_, sample_name]
    x <- rowSums(isd$insertSizes$genome[, .SD, .SDcols = stage_samples])
    x <- x / sum(x)    # express as a proportion of the total
    x <- pmax(1e-5, x) # prevent log(0) and impossible values
    log(x / sum(x))    # normalize to sum to 1 and take the log for NRLL calculation
}
extractInsertSizeEps <- function(isd, env){
    eps <- data.table(
        nucleosomal    = getStateEmissProbs(isd, env$NUCLEOSOMAL_STAGE),
        subnucleosomal = getStateEmissProbs(isd, env$SUBNUCLEOSOMAL_STAGE)
    )
    emissProbsFile <- paste(env$SHM_FILE_PREFIX, "emissionProbs_insertSize.tsv", sep = '.')
    write.table(
        eps, 
        file = emissProbsFile, 
        quote = FALSE, 
        row.names = FALSE, 
        col.names = FALSE, 
        sep = "\t"
    )
    emissProbsFile
}

# analyze and aggregate distributions of different bin scores
# all scores are expected to be one a comparable scale between samples, including
analyzeScoreDist <- function(bins, gcLimits, scores, distUnit = 0.1){
    scores_wrk <- scores[getIncudedAutosomeBins(bins, gcLimits)] # only use good autosomal bins to analyze score distributions
    dist <- data.table(x = as.integer(round(scores_wrk / distUnit)) * distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / replaceNaN(sum(y, na.rm = TRUE))]
    mu <- replaceNaN(mean(scores_wrk, na.rm = TRUE))
    sd <- replaceNaN(sd(  scores_wrk, na.rm = TRUE))
    list(
        score    = scores,
        dist     = dist,
        mean     = mu,
        sd       = sd,
        z        = (scores - mu) / sd,      # for normal/parametric scores
        quantile = ecdf(scores_wrk)(scores) # for non-parametric scores
    )
}
analyzeSampleScores <- function(bd, gcLimits, scoreFn, env, distUnit = 0.1, ...){
    x <- mclapply(bd$samples$sample_name, function(sample_name){
    # x <- lapply(bd$samples[filename_prefix %in% c("24290X11", "24290X9"), sample_name], function(sample_name){
        message(paste("   ", "analyzeSampleScores", sample_name))
        analyzeScoreDist(bd$bins$genome, gcLimits, scoreFn(bd, sample_name, ...), distUnit)
    }, mc.cores = env$N_CPU)
    # })
    names(x) <- bd$samples$sample_name
    x
}
aggregateAndAnalyzeScores <- function(bins, gcLimits, sampleScores, sample_names, distUnit = 0.1){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        sampleScores[[sample_name]]$score
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(bins, gcLimits, x[, replaceNaN(rowMeans(.SD, na.rm = TRUE))], distUnit)
}
aggregateSampleScores <- function(bd, gcLimits, stageTypes, sampleScores, env, distUnit = 0.1){

    # aggregate scores by spermatid stage
    allStages <- unique(bd$samples$stage)
    by_stage <- mclapply(allStages, function(stage_){
        message(paste("   ", "aggregateSampleScores by_stage", stage_))
        sample_names <- bd$samples[stage == stage_, sample_name]
        aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, distUnit)
    }, mc.cores = env$N_CPU)
    names(by_stage) <- allStages

    # aggregate scores by spermatid stage type (round vs. elong)
    by_stageType <- mclapply(names(stageTypes), function(stageType){
        message(paste("   ", "aggregateSampleScores by_stageType", stageType))
        sample_names <- bd$samples[stage %in% stageTypes[[stageType]], sample_name]
        aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, distUnit)
    }, mc.cores = env$N_CPU)
    names(by_stageType) <- names(stageTypes)

    # calculate the difference in scores between round and elongated spermatids to assess spermiogenesis trajectory
    stageType1 <- names(stageTypes)[1]
    stageType2 <- names(stageTypes)[2]
    stageType_delta <- analyzeScoreDist(
        bd$bins$genome, 
        gcLimits, 
        by_stageType[[stageType1]]$score - by_stageType[[stageType2]]$score, # typically round - elong to give positive values in early stages
        distUnit
    )
    list(
        by_stage        = by_stage,
        by_stageType    = by_stageType,
        stageType_delta = stageType_delta
    )
}

# # analyze and aggregate distributions of different bin scores
# # all scores are expected to be one a comparable scale between samples, including
# analyzeScoreDist <- function(bins, gcLimits, scores, distUnit = 0.1){
#     scores_wrk <- scores[getIncudedAutosomeBins(bins, gcLimits)] # only use good autosomal bins to analyze score distributions
#     dist <- data.table(x = as.integer(round(scores_wrk / distUnit)) * distUnit)[, .(y = .N), keyby = .(x)]
#     dist[, y := y / sum(y, na.rm = TRUE)]
#     mu <- mean(scores_wrk, na.rm = TRUE)
#     sd <- sd(  scores_wrk, na.rm = TRUE)
#     list(
#         score    = scores,
#         dist     = dist,
#         mean     = mu,
#         sd       = sd,
#         z        = (scores - mu) / sd,      # for normal/parametric scores
#         quantile = ecdf(scores_wrk)(scores) # for non-parametric scores
#     )
# }
# analyzeSampleScores <- function(bd, gcLimits, scoreFn, distUnit = 0.1){
#     # mclapply
#     sapply(bd$samples$sample_name, function(sample_name){
#         message(paste("   ", "analyzeSampleScores", sample_name))
#         analyzeScoreDist(bd$bins$genome, gcLimits, scoreFn(bd, sample_name), distUnit)
#     }, simplify = FALSE, USE.NAMES = TRUE)  
# }
# aggregateAndAnalyzeScores <- function(bins, gcLimits, sampleScores, sample_names, distUnit = 0.1){
#     # mclapply
#     x <- as.data.table(sapply(sample_names, function(sample_name){
#         sampleScores[[sample_name]]$score
#     }, simplify = FALSE, USE.NAMES = TRUE))
#     analyzeScoreDist(bins, gcLimits, x[, rowMeans(.SD, na.rm = TRUE)], distUnit) # TODO: this is the slowest step, can we speed it up?
# }
# aggregateSampleScores <- function(bd, gcLimits, stageTypes, sampleScores, distUnit = 0.1){
#     by_stage <- sapply(unique(bd$samples$stage), function(stage_){
#         sample_names <- bd$samples[stage == stage_, sample_name]
#         aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, distUnit)
#     }, simplify = FALSE, USE.NAMES = TRUE)
#     by_stageType <- sapply(names(stageTypes), function(stageType){
#         sample_names <- bd$samples[stage %in% stageTypes[[stageType]], sample_name]
#         aggregateAndAnalyzeScores(bd$bins$genome, gcLimits, sampleScores, sample_names, distUnit)
#     }, simplify = FALSE, USE.NAMES = TRUE)
#     stageType1 <- names(stageTypes)[1]
#     stageType2 <- names(stageTypes)[2]
#     stageType_delta <- analyzeScoreDist(
#         bd$bins$genome, 
#         gcLimits, 
#         by_stageType[[stageType1]]$score - by_stageType[[stageType2]]$score, # typically round - elong to give positive values in early stages
#         distUnit
#     )
#     list(
#         by_stage        = by_stage,
#         by_stageType    = by_stageType,
#         stageType_delta = stageType_delta
#     )
# }
