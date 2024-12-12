# establish the color palette for heat map and other visualizations
fold_change_color <- function(foldChange, maxFold){
    minFold <- 1 / maxFold 
    lfc <- log(pmax(minFold, pmin(maxFold, foldChange)), base = maxFold) # ranges from -1 to 1
    I <- floor(nTrackMapColorsPerSide * abs(lfc)) + 1L
    ifelse(lfc < 0, trackMapColors$low[I], trackMapColors$high[I])
}
z_score_color <- function(zScore, maxZScore){
    minZScore <- -maxZScore 
    z <- pmax(minZScore, pmin(maxZScore, zScore))
    I <- floor(nTrackMapColorsPerSide * abs(z) / maxZScore) + 1L
    ifelse(z < 0, trackMapColors$low[I], trackMapColors$high[I])
}
quantile_score_color <- function(quantile, maxQuantile){
    minQuantile <- 1 - maxQuantile
    quantile <- pmax(minQuantile, pmin(maxQuantile, quantile))
    I <- floor(nTrackMapColorsPerSide * abs(quantile - 0.5) / maxQuantile) + 1L
    ifelse(quantile < 0.5, trackMapColors$low[I], trackMapColors$high[I])
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
