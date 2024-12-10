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