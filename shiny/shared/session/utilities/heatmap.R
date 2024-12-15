# assemble the score heat map, e.g., for the track browser

# create a resuable horizontal black line as a row score space separator
scoreMapSeparatorImage <- function(config){
    x <- matrix("black", nrow = config$plotWidthPixels, ncol = config$sepHeightPixels) %>% 
    matrixToCImg(config$plotWidthPixels, config$sepHeightPixels)
    imager::imappend(
        list(
            scoreMapLabelImage(config$sepHeightPixels, config),
            x,
            scoreMapLegendImage(config$sepHeightPixels, config)
        ),
        axis = 'x'
    )
}

# create subimages as track image labels for the...
# ...left side label; this is generally the unit for the row display
scoreMapLabelImage <- function(height, config, label = NULL){
    x <- matrix("white", nrow = config$labelWidthPixels, ncol = height) %>% 
    matrixToCImg(config$labelWidthPixels, height)
    if(!is.null(label)) x <- x %>% imager::draw_text(5, 1, label, "black", fsize = height - 2)
    x
}
# ...right side label (in the legend area); this is generally the name of the sample or group
scoreMapLegendImage <- function(height, config, label = NULL){
    x <- matrix("white", nrow = config$legendWidthPixels, ncol = height) %>% 
    matrixToCImg(config$legendWidthPixels, height)
    if(!is.null(label)) x <- x %>% imager::draw_text(5, 1, label, "black", fsize = height - 2)
    x
}
# ...group header, representing a single score type; text identifies the score type and input data
scoreMapHeaderImage <- function(scoreType, config){
    x <- matrix("white", nrow = config$plotWidthPixels, ncol = config$headerHeightPixels) %>% 
    matrixToCImg(config$plotWidthPixels, config$headerHeightPixels) %>% 
    imager::draw_text(1, 2, scoreType$trackHeaderLabel, "black", fsize = config$headerHeightPixels - 2)
    paScoreMapYBreaks <<- rbind(paScoreMapYBreaks, data.table(
        height = config$headerHeightPixels + 1, scoreTypeName = "NA", rowType = "header", seriesName = "NA"
    ))
    imager::imappend(
        list(
            scoreMapLabelImage(config$headerHeightPixels, config),
            x,
            scoreMapLegendImage(config$headerHeightPixels, config)
        ),
        axis = 'x'
    )
}

# create subimages as track image rows representing one sample or group score value across display bins
scoreMapRowImage <- function(scoreTypeName, binI, b, scoreValues, colorFn, config, rowType, labelLabel = NULL, legendLabel = NULL){
    x <- merge(
        b, 
        data.table(
            binI = binI, 
            val  = scoreValues
        ),
        by = "binI", 
        all.x = TRUE
    )[, 
        .(val = weighted.mean(val, basesInPixel, na.rm = TRUE)), 
        keyby = .(pixel)
    ]
    if(nrow(x) != config$plotWidthPixels) x <- merge( # if the region extends beyond the chromosome boundary
        data.table(pixel = 1:config$plotWidthPixels),
        x,
        by = "pixel",
        all.x = TRUE
    )
    x <- colorFn(scoreTypeName, x$val, config) %>% 
    rep(config$Row_Height_Pixels) %>%
    matrixToCImg(config$plotWidthPixels, config$Row_Height_Pixels)
    seriesName <- if(is.null(legendLabel)) "NA" else legendLabel
    paScoreMapYBreaks <<- rbind(paScoreMapYBreaks, data.table(
        height = config$Row_Height_Pixels + 1, scoreTypeName = scoreTypeName, rowType = rowType, seriesName = seriesName
    ))
    imager::imappend(
        list(
            scoreMapLabelImage(config$Row_Height_Pixels, config, labelLabel),
            x,
            scoreMapLegendImage(config$Row_Height_Pixels, config, legendLabel)
        ),
        axis = 'x'
    )
}

# coordinate the assembly of a one complete score map group image for a single score type
# called by the paScoreMap track builder function
scoreMapGroupImage <- function(scoreTypeName, sourceId, bd, binI, b, config){
    startSpinner(session, message = paste("rendering", scoreTypeName))
    scoreLevel <- getScoreLevel(scoreTypeName)
    scoreType <- getScoreType(sourceId, scoreTypeName)
    sepImage <- scoreMapSeparatorImage(config)
    summaryScores <- getStageTypeDeltaScores(sourceId, scoreTypeName)[[1]][[scoreType$summaryType]][binI]
    seriesAggScores <- getSeriesAggScores(sourceId, scoreTypeName, bd$samples, config)
    imager::imappend(c(

        # a single row representing the fully aggregated score results for a scoreType
        # this is the only row for genome gc and txn scoreTypes
        list(
            scoreMapHeaderImage(scoreType, config),
            sepImage,
            scoreMapRowImage(
                scoreTypeName, binI, b, 
                summaryScores, 
                getSeriesSummaryColors, config, "summary", 
                labelLabel  = if(scoreLevel == "sample") "delta" else scoreType$trackSummaryLabel,
                legendLabel = if(scoreLevel == "sample") "round - elong" else NULL
            ),
            sepImage
        ),
        if(scoreLevel == "sample" && config$Aggregate_By != "none") c(

            # one or more rows representing sample-level primary scores
            # depending on the user setting for Aggregate_By, there may be one row per sample, stage, or stage type
            # these are _absolute_ scores; thus, it is possible for a single sample to have asymmetric score that are ~all high or low
            lapply(
                names(seriesAggScores),
                function(seriesName) scoreMapRowImage(
                    scoreTypeName, binI, b, 
                    seriesAggScores[[seriesName]]$score[binI], 
                    getSeriesSampleColors, config, "score",
                    labelLabel = if(names(seriesAggScores)[1] == seriesName) scoreType$trackScoreLabel else NULL,
                    legendLabel = seriesName
                )
            ),
            list(sepImage),

            # one or more rows representing intra-sample (or intra-group) _relative_ scores, expressed as per-sample/group quantiles
            # thus, these should show symmetric distributions for each sample or group with ~equal numbers of high and low scores
            # these rows are only not needed for GC residual Z scores which are already centered and symmetric
            if(scoreTypeName != "gcrz") lapply(
                names(seriesAggScores),
                function(seriesName) scoreMapRowImage(
                    scoreTypeName, binI, b, 
                    seriesAggScores[[seriesName]]$quantile[binI], 
                    getSeriesQuantileColors, config, "quantile",
                    labelLabel = if(names(seriesAggScores)[1] == seriesName) scoreType$summaryType else NULL,
                    legendLabel = seriesName
                )
            ) else list(),
            list(sepImage)
        ) else list()
    ), axis = 'y')
}
