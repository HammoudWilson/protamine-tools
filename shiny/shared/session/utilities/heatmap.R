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
scoreMapGroupImage <- function(scoreTypeName, sourceId, metadata, config, coord, binI, b){
    startSpinner(session, message = paste("rendering", scoreTypeName))
    scoreType <- getScoreType(scoreTypeName)
    scoreLevel <- getScoreLevel(scoreTypeName)
    isSampleScore <- scoreLevel == "sample"
    sepImage <- scoreMapSeparatorImage(config)
    scores <- list()
    summaryScores <- if(isSampleScore) {
        scores$quantile <- getSampleScores(metadata, config, coord, scoreTypeName, "quantile")
        scores$quantile$stageType # yes, this is stageType_delta (was named this way in score_functions.R)
    } else {
        getGenomeScores(sourceId, scoreTypeName, binI)
    }
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
        if(isSampleScore && config$Aggregate_By != "none") {
            seriesAggNames <- getSeriesAggNames(metadata, config)
            isGcrz <- startsWith(scoreTypeName, "gcrz")
            if(isGcrz && config$GCRZ_As_Quantiles) { # give user the option to plot GCRZ scores as quantiles
                scores$data <- getSampleScores(metadata, config, coord, scoreTypeName, "quantile")
                colorFn <- getSeriesQuantileColors
                labelLabel <- "quantile"
            } else {
                scores$data <- getSampleScores(metadata, config, coord, scoreTypeName, "score")
                colorFn <- getSeriesSampleColors
                labelLabel <- scoreType$trackScoreLabel
            }
            c(
                # one or more rows representing sample-level primary scores
                # depending on the user setting for Aggregate_By, there may be one row per sample, stage, or stage type
                # these may be _absolute_ scores; thus, a single sample may have asymmetric scores that are ~all high or low
                lapply(
                    seriesAggNames,
                    function(seriesName) scoreMapRowImage(
                        scoreTypeName, binI, b, 
                        scores$data[[seriesName]], 
                        colorFn, config, "score",
                        labelLabel = if(seriesAggNames[1] == seriesName) labelLabel else NULL,
                        legendLabel = seriesName
                    )
                ),
                list(sepImage),

                # one or more rows representing intra-sample (or intra-group) _relative_ scores, expressed as per-sample/group quantiles
                # thus, these should show symmetric distributions for each sample or group with ~equal numbers of high and low scores
                # these rows are not needed for GC residual Z scores which are already centered and symmetric relative values
                if(!isGcrz) {
                    if(is.null(scores$quantile)) scores$quantile <- getSampleScores(metadata, config, coord, scoreTypeName, "quantile")
                    lapply(
                        seriesAggNames,
                        function(seriesName) scoreMapRowImage(
                            scoreTypeName, binI, b,  
                            scores$quantile[[seriesName]], 
                            getSeriesQuantileColors, config, "quantile",
                            labelLabel = if(seriesAggNames[1] == seriesName) "quantile" else NULL,
                            legendLabel = seriesName
                        )
                    )
                } else list(),
                list(sepImage)
            ) 
        } else list()
    ), axis = 'y')
}
