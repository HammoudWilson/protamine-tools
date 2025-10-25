# utilities for parsing bin scores and distributions

# score types metadata; minimal information only as required to calculate score distributions
scoreTypes <- list(
    genome = list(
        gc = list(
            name = "gc",
            distUnit = 0.005,
            include = c("z"), # percent GC itself already available in bins object; this z harvested as gc_z
            log10 = FALSE
        ),
        txn = list(
            name = "txn",
            distUnit = 0.1,
            include = character(), # nascent transcription rate always handled on an absolute scale in app
            log10 = TRUE,
            minValue = 1e-3
        ),
        stgm = list(
            name = "stgm",
            distUnit = 0.1,
            include = c("quantile"), 
            log10 = FALSE
        ),
        hic = list(
            name = "hic",
            distUnit = 0.1 / 20,
            include = c("quantile"), 
            log10 = FALSE
        )
    ),
    sample = list(
        gcrz_obs = list(
            name = "gcrz_obs",
            # gcBiasDependent = TRUE, # thus, cannot be assessed until GC bias is established in app
            distUnit = 0.1,
            include = c("quantile"), # retain the ability to plot gcrz as values or relative quantiles
            log10 = FALSE            # they are supposed to be Z scores already but are not completely stable
        ),
        gcrz_wgt = list(
            name = "gcrz_wgt",
            # gcBiasDependent = TRUE,
            distUnit = 0.1,
            include = c("quantile"),
            log10 = FALSE
        ),
        nrll = list(
            name = "nrll",
            # gcBiasDependent = FALSE,
            distUnit = 0.1,
            include = c("quantile"), # retain the ability to plot NRLL as values or relative quantiles
            log10 = FALSE            # but these are mainly plotted on a custom scale in app
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

# extract the histone- and protamine-associated insert size distributions for NRLL calculation
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

# apply mappability and sometimes Tn5 site preference normalization to bin counts
normalize_mappability_tn5 <- function(sample_name_, normalizeTn5Site){

    # collect the counts or weights per bin per sample, based on normalizeTn5Site
    rpb <- if(normalizeTn5Site) collate$n_wgt_bin_smp[, sample_name_] 
                           else collate$n_obs_bin_smp[, sample_name_]

    # increase the counts of poorly mappable bins
    mpp_bin <- collate$mpp_bin_smp[, sample_name_]
    rpb <- ifelse(rpb > 0 & mpp_bin > 0, rpb / mpp_bin, 0)

    # rescale normalized counts to sum to the number of observed reads
    # this is done to maintain the same statistical weight pre- and post-normalization
    n_obs  <- sum(collate$n_obs_bin_smp[includedAutosomeBins, sample_name_])
    n_norm <- sum(rpb[includedAutosomeBins])
    rpb[isPrimaryGenome] <- rpb[isPrimaryGenome] / n_norm * n_obs
    rpb
}

# sample-level score type functions
# must return a vector the length of all composite bins, but scores only needed for primary genome bins
doGcRegression <- function(sample_name_, normalizeTn5Site){

    # collect bin data for this sample
    dt <- data.table(
        gc_bin = collate$gc_bin_smp[, sample_name_],
        rpb    = normalize_mappability_tn5(sample_name_, normalizeTn5Site) # reads per bin
    )

    # perform an initial GC regression to find a rough fit
    dt_fit <- dt[includedAutosomeBins, 
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
    dt_fit <- dt[includedAutosomeBins & zScore.nbinomCountsGC(fit, dt$rpb, dt$gc_bin, binCN = 2) %between% c(-3, 3)]
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
get_gcrz_obs <- function(sample_name_, ...){
    doGcRegression(sample_name_, normalizeTn5Site = FALSE)
}
get_gcrz_wgt <- function(sample_name_, ...){
    doGcRegression(sample_name_, normalizeTn5Site = TRUE)
}
get_nrll <- function(sample_name_, emissProbsFile){
    smp <- collate$samples[sample_name == sample_name_]
    script <- paste('bash', file.path(env$MODULES_DIR, 'score', 'get_bin_NRLL.sh'))
    fread(cmd = paste(script, smp$filename_prefix, emissProbsFile))[[1]] # one value per composite bin
}

# analyze and aggregate distributions of different bin scores
# all scores are expected to be one a comparable scale between samples
analyzeScoreDist <- function(scores, scoreType, data_set_name, return_scores = FALSE){

    # take log as needed, preventing log(0) and impossible values
    if(scoreType$log10) scores <- log10(pmax(scoreType$minValue, scores)) 
    
    # use only good autosomal bins to analyze score distributions
    scores_wrk <- scores[includedAutosomeBins] 
    dist <- data.table(x = as.integer(round(scores_wrk / scoreType$distUnit)) * scoreType$distUnit)[, .(y = .N), keyby = .(x)]
    dist[, y := y / replaceNaN(sum(y, na.rm = TRUE))]
    median <- replaceNaN(median(scores_wrk, na.rm = TRUE))
    mean   <- replaceNaN(mean(  scores_wrk, na.rm = TRUE))
    sd     <- replaceNaN(sd(    scores_wrk, na.rm = TRUE))

    # write a data.table to disk with the different bin-level scores
    scoreFile <- paste(env$TMP_FILE_PREFIX, "bin_scores", scoreType$name, data_set_name, "rds", sep = '.')
    saveRDS(data.table(
        score    = scores,
        z        = if("z"        %in% scoreType$include) (scores - median) / sd   else NULL, # for normal/parametric scores
        quantile = if("quantile" %in% scoreType$include) ecdf(scores_wrk)(scores) else NULL  # for non-parametric scores
    )[isPrimaryGenome], file = scoreFile)

    # return a list with the distribution and summary statistics
    list(
        dist   = dist,
        peak   = dist$x[which.max(dist$y)],
        median = median,
        mean   = mean,
        sd     = sd,
        scoreFile = scoreFile,
        scores = if(return_scores) scores else NULL # temporarily retain sample level scores for later aggregation
    )
}
analyzeSampleScores <- function(scoreFn, scoreType, ...){
    x <- mclapply(collate$samples$sample_name, function(sample_name){
    # x <- lapply(collate$samples[filename_prefix %in% c("24290X11", "24290X9"), sample_name], function(sample_name){
        message(paste("   ", "analyzeSampleScores", sample_name))
        analyzeScoreDist(scoreFn(sample_name, ...), scoreType, sample_name, return_scores = TRUE)
    }, mc.cores = env$N_CPU)
    # })
    names(x) <- collate$samples$sample_name
    x
}
aggregateAndAnalyzeScores <- function(sampleScores, sample_names, scoreType, data_set_name, return_scores = FALSE){
    x <- as.data.table(sapply(sample_names, function(sample_name){
        sampleScores[[sample_name]]$scores
    }, simplify = FALSE, USE.NAMES = TRUE))
    analyzeScoreDist(x[, replaceNaN(rowMeans(.SD, na.rm = TRUE))], scoreType, data_set_name, return_scores)
}
aggregateSampleScores <- function(sampleScores, scoreType){

    # aggregate scores by spermatid stage
    allStages <- unique(collate$samples$stage)
    by_stage <- mclapply(allStages, function(stage_){
    # by_stage <- lapply(allStages, function(stage_){
        message(paste("   ", "aggregateSampleScores by_stage", stage_))
        sample_names <- collate$samples[stage == stage_, sample_name]
        aggregateAndAnalyzeScores(sampleScores, sample_names, scoreType, stage_)
    }, mc.cores = env$N_CPU)
    # })
    names(by_stage) <- allStages

    # aggregate scores by spermatid stage type (round vs. elong)
    by_stageType <- mclapply(names(stageTypes), function(stageType){
    # by_stageType <- lapply(names(stageTypes), function(stageType){
        message(paste("   ", "aggregateSampleScores by_stageType", stageType))
        sample_names <- collate$samples[stage %in% stageTypes[[stageType]], sample_name]
        aggregateAndAnalyzeScores(
            sampleScores, sample_names, scoreType, stageType, 
            return_scores = stageType %in% names(stageTypes)[1:2]
        )
    }, mc.cores = env$N_CPU)
    # })
    names(by_stageType) <- names(stageTypes)

    # calculate the difference in scores between round and elongated spermatids to assess spermiogenesis trajectory
    message(paste("   ", "aggregateSampleScores stageType_delta"))
    stageType1 <- names(stageTypes)[1]
    stageType2 <- names(stageTypes)[2]
    delta <- list(
        stageType = analyzeScoreDist(
            by_stageType[[stageType1]]$scores - by_stageType[[stageType2]]$scores, # typically round - elong to give positive values in early stages
            scoreType,
            "stageType_delta"
        )
    )
    by_stageType[[stageType1]]$scores <- NULL
    by_stageType[[stageType2]]$scores <- NULL
    list(
        by_stage     = by_stage,
        by_stageType = by_stageType,
        delta        = delta
    )
}
