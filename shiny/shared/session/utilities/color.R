# establish the dynamic color palettes for heat map and other visualizations

# functions that convert different types of normalized scores to color ranges
z_score_color <- function(zScore, maxZScore){
    minZScore <- -maxZScore 
    z <- pmax(minZScore, pmin(maxZScore, zScore))
    I <- floor(nTrackMapColorsPerSide * abs(z) / maxZScore) + 1L
    ifelse(z < 0, trackMapColors$low[I], trackMapColors$high[I])
}
quantile_score_color <- function(quantile, maxQuantile){
    minQuantile <- 1 - maxQuantile
    quantile <- pmax(minQuantile, pmin(maxQuantile, quantile))
    I <- floor(nTrackMapColorsPerSide * abs(quantile - 0.5) / (maxQuantile - 0.5)) + 1L
    ifelse(quantile < 0.5, trackMapColors$low[I], trackMapColors$high[I])
}
cpm_score_color <- function(log10cpm, minLog10Cpm = -1, maxLog10Cpm = 1){
    log10cpm <- pmax(minLog10Cpm, pmin(maxLog10Cpm, log10cpm))
    I <- floor(nTrackMapColorsPerSide * (log10cpm - minLog10Cpm) / (maxLog10Cpm - minLog10Cpm)) + 1L
    trackMapColors$high[I]
}
nrll_score_color <- function(nrll, maxNrll){
    minNrll <- -maxNrll
    nrll <- pmax(minNrll, pmin(maxNrll, nrll))
    I <- floor(nTrackMapColorsPerSide * abs(nrll) / maxNrll) + 1L
    ifelse(nrll < 0, trackMapColors$low[I], trackMapColors$high[I])
}
fraction_score_color <- function(fraction, minFraction = 0, maxFraction = 1){
    midpoint <- (minFraction + maxFraction) / 2
    fraction <- pmax(minFraction, pmin(maxFraction, fraction))
    I <- floor(nTrackMapColorsPerSide * abs(fraction - midpoint) / (maxFraction - midpoint)) + 1L
    ifelse(fraction < midpoint, trackMapColors$low[I], trackMapColors$high[I])
}

# functions to get distribution trace colors
getSampleColorsByStage <- function(allSamples, samples){
    stages <- allSamples[, unique(stage)]
    colors <- sapply(allSamples$stage, function(x) stageColors[which(stages == x)])
    names(colors) <- allSamples$sample_name
    colors[samples$sample_name]
}
getStageColors <- function(allSamples, samples){
    stages <- allSamples[, unique(stage)]
    colors <- stageColors[1:length(stages)]
    names(colors) <- stages
    colors[samples[, unique(stage)]]
}
getStageTypeColors <- function(sourceId, allSamples, samples){
    allStageTypes <- unique(getStageTypesByStage(sourceId, allSamples[, unique(stage)]))
    stageTypes    <- unique(getStageTypesByStage(sourceId,    samples[, unique(stage)]))
    colors <- stageTypeColors[1:length(allStageTypes)]
    names(colors) <- allStageTypes
    colors[stageTypes]
}
getSeriesSummaryColors <- function(scoreTypeName, scoreValues, config){ # used to color the top group-level summary row of every track heatmap group
    switch(
        scoreTypeName,
        gc   = z_score_color(       scoreValues, config$Max_Z_Score),
        txn  = cpm_score_color(     scoreValues, config$Min_Txn_Log10_CPM, config$Max_Txn_Log10_CPM),
        gcrz = quantile_score_color(scoreValues, config$Max_Quantile), # all sample-level summary scores use quantiles, i.e., assume non-parametric distributions
        iisf = quantile_score_color(scoreValues, config$Max_Quantile),
        nrll = quantile_score_color(scoreValues, config$Max_Quantile)
    )
}
getSeriesSampleColors <- function(scoreTypeName, scoreValues, config){ # used to color the subsequent sample-level rows of every track heatmap group
    switch(
        scoreTypeName,
        gcrz = z_score_color(scoreValues,           config$Max_Z_Score),
        iisf = fraction_score_color(scoreValues, 0, config$Max_Fraction_IIS), # unlike gcrz and nrll, iisf is not inherently centered or symmetric
        nrll = nrll_score_color(scoreValues,        config$Max_NRLL)
    )
}
getSeriesQuantileColors <- function(scoreTypeName, scoreValues, config){ # like above, now coloring the quantile rows for intra-sample/group relative scores
    quantile_score_color(scoreValues, config$Max_Quantile)
}
