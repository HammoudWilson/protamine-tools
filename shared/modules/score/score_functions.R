# utilities for parsing bin scores and distributions

# score types metadata; minimal information only as required to calculate score distributions
scoreTypes <- list(
    genome = list(
        gc = list(
            distUnit = 0.005,
            include = c("z"), # percent GC itself already available in bins object
            log10 = FALSE
        ),
        txn = list(
            distUnit = 0.1,
            include = c("score"), # nascent transcription rate always handled on an absolute scale
            log10 = TRUE,
            minValue = 1e-3
        )
    ),
    sample = list(
        gcrz_obs = list(
            # gcBiasDependent = TRUE, # thus, cannot be assessed until GC bias is established in app
            distUnit = 0.1,
            include = c("score","z","quantile"), # z/quantile only relevant for delta (others already a Z)
            log10 = FALSE
        ),
        gcrz_wgt = list(
            # gcBiasDependent = TRUE,
            distUnit = 0.1,
            include = c("score","z","quantile"), # z/quantile only relevant for delta (others already a Z)
            log10 = FALSE
        ),
        nrll = list(
            # gcBiasDependent = FALSE,
            distUnit = 0.1,
            include = c("score","quantile"), # retain the ability to plot NRLL as values or relative quantiles
            log10 = FALSE
        )
    )
)

# score utility functions
unpackStageTypes <- function(env){ # unpack the stage types from the user option value
    stageTypes <- list()
    for(stageType in strsplit(env$STAGE_TYPES, ';')[[1]]){
        x <- strsplit(stageType, ':')[[1]]
        stageTypes[[x[1]]] <- strsplit(x[2], ',')[[1]]
    }
    stageTypes
}
replaceNaN <- function(x){
    x[is.nan(x)] <- NA
    x
}

# apply mappability and sometimes Tn5 site preference normalization
normalize_mappability_tn5 <- function(collate, sample_name_, gcLimits, normalizeTn5Site){

    # collect the counts or weights per bin per sample, based on normalizeTn5Site
    rpb <- if(normalizeTn5Site) collate$n_wgt_bin_smp[, sample_name_] 
                           else collate$n_obs_bin_smp[, sample_name_]

    # increase the counts of poorly mappable bins
    mpp_bin <- collate$mpp_bin_smp[, sample_name_]
    rpb <- ifelse(rpb > 0 & mpp_bin > 0, rpb / mpp_bin, 0)

    # rescale normalized counts to sum to the number of observed reads
    # this is done to maintain the same statistical weight pre- and post-normalization
    for (genome in c(collate$env$PRIMARY_GENOME, collate$env$SPIKE_IN_GENOME)) {
        I_ref  <- getGenomeBins(collate$bins, genome)
        I_norm <- getIncludedAutosomeBins(collate$bins, genome, gcLimits)
        n_obs  <- sum(collate$n_obs_bin_smp[I_norm, sample_name_])
        n_norm <- sum(rpb[I_norm])
        rpb[I_ref] <- rpb[I_ref] / n_norm * n_obs
    }
    rpb
}

# sample-level score type functions
doGcRegression <- function(collate, sample_name_, gcLimits, normalizeTn5Site){

    # determine the set of included autosomal bins used for the regression
    fitI <- getIncludedAutosomeBins(collate$bins, collate$env$PRIMARY_GENOME, gcLimits)

    # collect bin data for this sample
    dt <- data.table(
        gc_bin = collate$gc_bin_smp[, sample_name_],
        rpb    = normalize_mappability_tn5(collate, sample_name_, gcLimits, normalizeTn5Site) # reads per bin
    )

    # perform an initial GC regression to find a rough fit
    dt_fit <- dt[fitI, 
        .(
            q025 = quantile(rpb, 0.025), # use a rough 2.5% to 97.5% quantile approximation to exclude outliers intially
            q975 = quantile(rpb, 0.975),
            gc_bin,
            rpb
        ), 
        keyby = .(x = as.integer(round(gc_bin * 20)))
    ][rpb >= q025 & rpb <= q975]
    J <- sample(1:nrow(dt_fit), 200000)
    fit <- new_nbinomCountsGC(
        binCounts  = dt_fit$rpb[J], # downsample the rough fit for a bit of speed
        fractionGC = dt_fit$gc_bin[J],
        binCN      = 2, # since working on autosomal bins per above
        method     = 'cubic'
    )

    # now use the rough fit to remove outliers > 3 std. dev. and re-fit
    dt_fit <- dt[fitI & zScore.nbinomCountsGC(fit, dt$rpb, dt$gc_bin, binCN = 2) %between% c(-3, 3)]
    fit <- new_nbinomCountsGC(
        binCounts  = dt_fit$rpb, # this fit uses all included data points
        fractionGC = dt_fit$gc_bin,
        binCN      = 2,
        method     = 'cubic'
    )

    # save this fit to a tmp file for later assembly into a list of all samples
    grcz_type <- if(normalizeTn5Site) "gcrz_wgt" else "gcrz_obs"
    tmpFile <- paste(env$SHM_FILE_PREFIX, "gcBiasModel", grcz_type, sample_name_, "rds", sep = '.')
    saveRDS(fit, file = tmpFile)

    # calculate the final z-scores for the sample as the gcrz_xxx score type
    z <- zScore.nbinomCountsGC(fit, dt$rpb, dt$gc_bin, binCN = collate$bins$nAlleles)
    z[abs(z) == Inf] <- NA
    z
}
get_gcrz_obs <- function(collate, sample_name_, gcLimits, ...){
    doGcRegression(collate, sample_name_, gcLimits, normalizeTn5Site = FALSE)
}
gcrz_wgt <- function(collate, sample_name_, gcLimits, ...){
    doGcRegression(collate, sample_name_, gcLimits, normalizeTn5Site = TRUE)
}
get_nrll <- function(collate, sample_name_, gcLimits, emissProbsFile){
    smp <- collate$samples[sample_name == sample_name_]
    script <- paste('bash', file.path(env$ACTION_DIR, 'get_bin_NRLL.sh'))
    fread(cmd = paste(script, smp$filename_prefix, emissProbsFile))[[1]] # one value per composite bin
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
    scores_wrk <- scores[getIncludedAutosomeBins(bins, genome, gcLimits)] # only use good autosomal bins to analyze score distributions
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
        analyzeScoreDist(collate$bins, genome, gcLimits, scoreFn(collate, sample_name, gcLimits, ...), scoreType)
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
