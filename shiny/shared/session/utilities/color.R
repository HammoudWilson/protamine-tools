# establish the color palettes for heat map and other visualizations

#----------------------------------------------------------------------
# get (staged) distribution trace colors for consistent plotting
#----------------------------------------------------------------------
stageColors <- c(
    CONSTANTS$plotlyColors$black, # spermatid stages, e.g., early_RS, etc.
    CONSTANTS$plotlyColors$blue,
    CONSTANTS$plotlyColors$orange,
    CONSTANTS$plotlyColors$green,
    CONSTANTS$plotlyColors$purple,
    CONSTANTS$plotlyColors$red,
    CONSTANTS$plotlyColors$yellow,
    CONSTANTS$plotlyColors$teal
)
stageTypeColors <- c(
    CONSTANTS$plotlyColors$red, # spermatid stage types, e.g., round, etc.
    CONSTANTS$plotlyColors$blue,
    CONSTANTS$plotlyColors$purple
)
getSampleColorsByStage <- function(allSamples, samples){
    orderedSamples <- allSamples[order(staging_order)]
    stages <- orderedSamples[, unique(stage)]
    colors <- sapply(orderedSamples$stage, function(x) stageColors[which(stages == x)])
    names(colors) <- orderedSamples$sample_name
    colors[samples$sample_name]
}
getSampleColorsBySample <- function(samples){
    orderedSamples <- samples[order(staging_order)]
    colors <- stageColors[1:nrow(orderedSamples)]
    names(colors) <- orderedSamples$sample_name
    colors
}
getStageColor <- function(metadata, stage){
    orderedSamples <- metadata$samples[order(staging_order)]
    stages <- orderedSamples[, unique(stage)]
    colors <- stageColors[1:length(stages)]
    names(colors) <- stages
    colors[stage]
}
getStageColors <- function(allSamples, samples){
    orderedSamples <- allSamples[order(staging_order)]
    stages <- orderedSamples[, unique(stage)]
    colors <- stageColors[1:length(stages)]
    names(colors) <- stages
    colors[samples[, unique(stage)]]
}
getStageTypeColor <- function(metadata, stageType){
    stageTypeStages <- metadata$stageTypes[[stageType]]
    getStageColor(metadata, stageTypeStages[1])
}
getStageTypeColors <- function(sourceId, allSamples, samples){
    orderedSamples <- allSamples[order(staging_order)]
    allStageTypes <- unique(getStageTypesByStage(sourceId, orderedSamples[, unique(stage)]))
    stageTypes    <- unique(getStageTypesByStage(sourceId,        samples[, unique(stage)]))
    colors <- stageTypeColors[1:length(allStageTypes)]
    names(colors) <- allStageTypes
    colors[stageTypes]
}

#----------------------------------------------------------------------
# convert different types of normalized scores to dynamic color ranges for heat maps
#----------------------------------------------------------------------

# establish a range of 61 colors for heat maps, with 30 colors on each side of the neutral color (grey)
paColors <- list(
    RED  = rgb(0.9, 0,   0),
    GREY = rgb(0.75, 0.75, 0.75),
    BLUE = rgb(0,   0,   1)
)
nTrackMapColorsPerSide <- 30
trackMapColors <- list(
    low  = colorRampPalette(c(paColors$GREY, paColors$BLUE))(nTrackMapColorsPerSide + 1), # blue color is cold/depleted,
    high = colorRampPalette(c(paColors$GREY, paColors$RED))( nTrackMapColorsPerSide + 1)  # red  color is hot/ enriched
)

# implement color spreading functions for different score types
z_score_color <- function(zScore, config){
    minZScore <- -config$Max_Z_Score
    z <- pmax(minZScore, pmin(config$Max_Z_Score, zScore))
    I <- floor(nTrackMapColorsPerSide * abs(z) / config$Max_Z_Score) + 1L
    ifelse(z < 0, trackMapColors$low[I], trackMapColors$high[I])
}
# quantile_score_color <- function(quantile, config){
#     minQuantile <- 1 - maxQuantile
#     quantile <- pmax(minQuantile, pmin(maxQuantile, quantile))
#     I <- floor(nTrackMapColorsPerSide * abs(quantile - 0.5) / (maxQuantile - 0.5)) + 1L
#     ifelse(quantile < 0.5, trackMapColors$low[I], trackMapColors$high[I])
# }
quantile_score_color <- function(quantile, config){
    # prevent infinite Z-scores by clamping extreme quantiles to ~6 sigma, i.e., 1e-9
    quantile <- pmax(1e-9, pmin(1-1e-9, quantile))
    # use normal distribution to weight quantiles
    z <- qnorm(quantile)
    z_score_color(z, config)
}
cpm_score_color <- function(log10cpm, config){
    min <- config$Min_Txn_Log10_CPM
    max <- config$Max_Txn_Log10_CPM
    log10cpm <- pmax(min, pmin(max, log10cpm))
    I <- floor(nTrackMapColorsPerSide * (log10cpm - min) / (max - min)) + 1L
    trackMapColors$high[I]
}
nrll_score_color <- function(nrll, config){
    NRLL_Limits <- as.numeric(strsplit(config$NRLL_Limits, ",")[[1]])
    low  <- NRLL_Limits[1]
    mid  <- NRLL_Limits[2]
    high <- NRLL_Limits[3]
    nrll <- pmax(low, pmin(high, nrll))
    I <- ifelse(nrll < mid,
        # scale in lower range
        floor(nTrackMapColorsPerSide * abs(nrll - mid) / abs(low - mid)) + 1L,
        # scale in upper range (enriched for small insert pattern)
        floor(nTrackMapColorsPerSide * abs(nrll - mid) / abs(high - mid)) + 1L
    )
    ifelse(nrll < mid, trackMapColors$low[I], trackMapColors$high[I])
}
fraction_score_color <- function(fraction, config){
    minFraction <- config$Min_Fraction
    maxFraction <- config$Max_Fraction
    midpoint <- (minFraction + maxFraction) / 2
    fraction <- pmax(minFraction, pmin(maxFraction, fraction))
    I <- floor(nTrackMapColorsPerSide * abs(fraction - midpoint) / (maxFraction - midpoint)) + 1L
    ifelse(fraction < midpoint, trackMapColors$low[I], trackMapColors$high[I])
}

#----------------------------------------------------------------------
# get heatmap track colors for different score types and aggregation levels 
#----------------------------------------------------------------------

# color the top group-level summary row of every track heatmap group
# this is the single score value for genome-level scores
# or the stageType delta for sample-level scores
getSeriesSummaryColors <- function(scoreTypeName, scoreValues, config){ # used to color the top group-level summary row of every track heatmap group
    switch(
        scoreTypeName,

        # genome scores here are equivalent to sample-level getSeriesSampleColors below
        # provided by scoreMapGroupImage as raw score values, used as is
        gc_z = z_score_color(  scoreValues, config), 
        txn  = cpm_score_color(scoreValues, config),

        # sample-level summary scores (deltas) use quantiles, i.e., assume non-parametric distributions
        # always provided as quantiles by scoreMapGroupImage
        quantile_score_color(scoreValues, config)
    )
}

# color the subsequent sample-level rows of every track heatmap group
# i.e., these values are provided as absolute scores by scoreMapGroupImage
getSeriesSampleColors <- function(scoreTypeName, scoreValues, config){ 
    switch(
        scoreTypeName,
        gcrz_obs = z_score_color(   scoreValues, config),  # these score values are already Z scores
        gcrz_wgt = z_score_color(   scoreValues, config),
        nrll     = nrll_score_color(scoreValues, config)
    )
}

# color final quantile rows for intra-sample/intra-group relative scores
# these values are always provided as quantiles by scoreMapGroupImage
# also used for GCRZ scores if the user requested quantiles
getSeriesQuantileColors <- function(scoreTypeName, scoreValues, config){  # scoreTypeName unused but needed for consistency
    quantile_score_color(scoreValues, config)
}
