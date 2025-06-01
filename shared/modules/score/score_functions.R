# utilities for parsing bin scores and distributions

# sample-level score type functions
replaceNaN <- function(x){
    x[is.nan(x)] <- NA
    x
}
get_nrll <- function(collate, sample_name_, emissProbsFile){
    smp <- collate$samples[sample_name == sample_name_]
    script <- paste('bash', file.path(env$ACTION_DIR, 'get_bin_NRLL.sh'))
    fread(cmd = paste(script, smp$filename_prefix, emissProbsFile))[[1]] # one value per composite bin
}

# score types metadata; minimal information only as required to calculate score distributions
scoreTypes <- list(
    genome = list(
        gc = list(
            gcBiasDependent = FALSE,
            distUnit = 0.005,
            include = c("z"), # percent GC itself already available in bins object
            log10 = FALSE
        ),
        txn = list(
            gcBiasDependent = FALSE,
            distUnit = 0.1,
            include = c("score"), # nascent transcription rate always handled on an absolute scale
            log10 = TRUE,
            minValue = 1e-3
        )
    ),
    sample = list(
        gcrz = list(
            gcBiasDependent = TRUE, # thus, cannot be assessed until GC bias is established in app
            distUnit = 0.1,
            include = c("score","z","quantile"), # z/quantile only relevant for delta (others already a Z); not calculated here in any case
            log10 = FALSE
        ),
        nrll = list(
            gcBiasDependent = FALSE,
            distUnit = 0.1,
            include = c("score","quantile"), # retain the ability to plot NRLL as values or relative quantiles
            log10 = FALSE
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

# extract the histone- and protamine-associated insert size distributions
# where a single specific stage is taken as being a sufficiently pure representation of a state
getStateEmissProbs <- function(samples, f_obs_isl_smp, stage_){ 
    stage_samples <- samples[stage == stage_, sample_name]
    f_obs_isl <- rowSums(f_obs_isl_smp[, stage_samples])
    f_obs_isl <- f_obs_isl / sum(f_obs_isl) # express as a proportion of the total
    f_obs_isl <- pmax(1e-5, f_obs_isl)      # prevent log(0) and impossible values
    log(f_obs_isl / sum(f_obs_isl))         # normalize to sum to 1 and take the log for NRLL calculation
}
extractInsertSizeEps <- function(samples, f_obs_isl_smp, env){
    eps <- data.table(
        histone   = getStateEmissProbs(samples, f_obs_isl_smp, env$HISTONE_STAGE),
        protamine = getStateEmissProbs(samples, f_obs_isl_smp, env$PROTAMINE_STAGE)
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
analyzeScoreDist <- function(bins, genome, gcLimits, scores, scoreType){
    if(scoreType$log10) scores <- log10(pmax(scoreType$minValue, scores)) # prevent log(0) and impossible values
    scores_wrk <- scores[getIncudedAutosomeBins(bins, genome, gcLimits)] # only use good autosomal bins to analyze score distributions
    dist <- data.table(x = as.integer(round(scores_wrk / scoreType$distUnit)) * scoreType$distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / replaceNaN(sum(y, na.rm = TRUE))]
    mu <- replaceNaN(mean(scores_wrk, na.rm = TRUE))
    sd <- replaceNaN(sd(  scores_wrk, na.rm = TRUE))
    median <- replaceNaN(median(scores_wrk, na.rm = TRUE))
    list(
        score    = scores, # score are dropped downstream when not included in scoreType
        dist     = dist,
        median   = median,
        mean     = mu,
        sd       = sd,
        peak     = dist$x[which.max(dist$y)],
        z        = if("z" %in% scoreType$include) (scores - median) / sd else NULL,         # for normal/parametric scores
        quantile = if("quantile" %in% scoreType$include) ecdf(scores_wrk)(scores) else NULL # for non-parametric scores
    )
}
analyzeSampleScores <- function(collate, genome, gcLimits, scoreFn, env, scoreType, ...){
    x <- mclapply(collate$samples$sample_name, function(sample_name){
    # x <- lapply(collate$samples[filename_prefix %in% c("24290X11", "24290X9"), sample_name], function(sample_name){
        message(paste("   ", "analyzeSampleScores", sample_name))
        analyzeScoreDist(collate$bins, genome, gcLimits, scoreFn(collate, sample_name, ...), scoreType)
    }, mc.cores = env$N_CPU)
    # })
    names(x) <- collate$samples$sample_name
    x
}
aggregateAndAnalyzeScores <- function(bins, genome, gcLimits, sampleScores, sample_names, scoreType){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        sampleScores[[sample_name]]$score
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(bins, genome, gcLimits, x[, replaceNaN(rowMeans(.SD, na.rm = TRUE))], scoreType)
}
aggregateSampleScores <- function(collate, genome, gcLimits, stageTypes, sampleScores, env, scoreType){

    # aggregate scores by spermatid stage
    allStages <- unique(collate$samples$stage)
    by_stage <- mclapply(allStages, function(stage_){
    # by_stage <- lapply(allStages, function(stage_){
        message(paste("   ", "aggregateSampleScores by_stage", stage_))
        sample_names <- collate$samples[stage == stage_, sample_name]
        aggregateAndAnalyzeScores(collate$bins, genome, gcLimits, sampleScores, sample_names, scoreType)
    }, mc.cores = env$N_CPU)
    # })
    names(by_stage) <- allStages

    # aggregate scores by spermatid stage type (round vs. elong)
    by_stageType <- mclapply(names(stageTypes), function(stageType){
    # by_stageType <- lapply(names(stageTypes), function(stageType){
        message(paste("   ", "aggregateSampleScores by_stageType", stageType))
        sample_names <- collate$samples[stage %in% stageTypes[[stageType]], sample_name]
        aggregateAndAnalyzeScores(collate$bins, genome, gcLimits, sampleScores, sample_names, scoreType)
    }, mc.cores = env$N_CPU)
    # })
    names(by_stageType) <- names(stageTypes)

    # calculate the difference in scores between round and elongated spermatids to assess spermiogenesis trajectory
    stageType1 <- names(stageTypes)[1]
    stageType2 <- names(stageTypes)[2]
    stageType_delta <- analyzeScoreDist(
        collate$bins,
        genome,  
        gcLimits, 
        by_stageType[[stageType1]]$score - by_stageType[[stageType2]]$score, # typically round - elong to give positive values in early stages
        scoreType
    )
    list(
        by_stage        = by_stage,
        by_stageType    = by_stageType,
        stageType_delta = stageType_delta
    )
}
