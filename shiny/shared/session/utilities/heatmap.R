# assemble the score heat map, e.g., for the track browser

# create a resuable horizontal black line as a row score space separator
scoreMapSeparatorImage <- function(config){
    x <- matrix("black", nrow = config$plotWidthPixels, ncol = config$sepHeightPixels) %>% 
    matrixToCImg(config$plotWidthPixels, config$sepHeightPixels)
    imager::imappend(
        list(
            scoreMapLabelImage(config$sepHeightPixels, config$sepHeightPixels, config),
            x,
            scoreMapLegendImage(config$sepHeightPixels, config$sepHeightPixels, config)
        ),
        axis = 'x'
    )
}

# create subimages as track image labels for the...
# ...left side label; this is generally the unit for the row display
scoreMapLabelImage <- function(trackHeight, labelHeight, config, label = NULL){
    x <- matrix("white", nrow = config$labelWidthPixels, ncol = trackHeight) %>% 
    matrixToCImg(config$labelWidthPixels, trackHeight)
    if(!is.null(label)) x <- x %>% imager::draw_text(5, 1, label, "black", fsize = labelHeight - 2)
    x
}
# ...right side label (in the legend area); this is generally the name of the sample or group
scoreMapLegendImage <- function(trackHeight, labelHeight, config, label = NULL){
    x <- matrix("white", nrow = config$legendWidthPixels, ncol = trackHeight) %>% 
    matrixToCImg(config$legendWidthPixels, trackHeight)
    if(!is.null(label)) x <- x %>% imager::draw_text(5, 1, label, "black", fsize = labelHeight - 2)
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
            scoreMapLabelImage(config$headerHeightPixels, config$headerHeightPixels, config),
            x,
            scoreMapLegendImage(config$headerHeightPixels, config$headerHeightPixels, config)
        ),
        axis = 'x'
    )
}

# create subimages as track image rows representing one sample or group score value across display bins
scoreMapRowImage <- function(scoreTypeName, isCutTagScore, isLog10, binI, b, scoreValues, 
                             colorFn, config, rowType, labelLabel = NULL, legendLabel = NULL){
    isCutTagHeatmap <- isCutTagScore && config$CutTag_As_HeatMap
    rowHeightPixels <- if(isCutTagScore){
        k_bin   <- config$binSize / 1e3
        m_reads <- config$samples[
            sample_name == legendLabel | stage == legendLabel, 
            sum(n_reads) / 1e6
        ]
        if(isCutTagHeatmap) config$Row_Height_Pixels else config$Row_Height_Pixels * 2
    } else {
        config$Row_Height_Pixels
    }
    x <- merge(
        b, 
        data.table(
            binI = binI, 
            val  = scoreValues
        ),
        by = "binI", 
        all.x = TRUE
    )[, 
        .(
            val = if(isCutTagScore)      weighted.mean(val / k_bin / m_reads, basesInPixel, na.rm = TRUE)
                  else if(isLog10) log10(weighted.mean(10**val,               basesInPixel, na.rm = TRUE))
                  else                   weighted.mean(    val,               basesInPixel, na.rm = TRUE)
        ), 
        keyby = .(pixel)
    ]
    if(nrow(x) != config$plotWidthPixels) x <- merge( # if the region extends beyond the chromosome boundary
        data.table(pixel = 1:config$plotWidthPixels),
        x,
        by = "pixel",
        all.x = TRUE
    )
    x <- if(isCutTagScore && !isCutTagHeatmap){
        nullBin <- rep("white", rowHeightPixels)
        ymax <- quantile(x$val, 0.975, na.rm = TRUE)
        y <- as.integer(round(pmin(ymax, x$val) / ymax * rowHeightPixels))
        as.vector(t(sapply(y, function(y) {
            if(is.na(y)) nullBin 
            else c(rep("white", rowHeightPixels - y), rep("black", y))
        })))
    } else {
        colorFn(scoreTypeName, x$val, config) %>% 
        rep(rowHeightPixels)
    }

    seriesName <- if(is.null(legendLabel)) "NA" else legendLabel
    paScoreMapYBreaks <<- rbind(paScoreMapYBreaks, data.table(
        height = rowHeightPixels + 1, scoreTypeName = scoreTypeName, rowType = rowType, seriesName = seriesName
    ))
    imager::imappend(
        list(
            scoreMapLabelImage(rowHeightPixels, config$Row_Height_Pixels, config, labelLabel),
            matrixToCImg(x, config$plotWidthPixels, rowHeightPixels),
            scoreMapLegendImage(rowHeightPixels, config$Row_Height_Pixels, config, legendLabel)
        ),
        axis = 'x'
    )
}

# coordinate the assembly of a one complete score map group image for a single score type
# called by the paScoreMap track builder function
scoreMapGroupImage <- function(scoreTypeName, sourceId, atac_metadata, cuttag_metadata, config, coord, binI, b){
    startSpinner(session, message = paste("rendering", scoreTypeName))
    scoreType <- getScoreType(scoreTypeName)
    scoreLevel <- getScoreLevel(scoreTypeName)
    isSampleScore <- scoreLevel == "sample"
    isCutTagScore <- !is.null(scoreType$cuttag) && scoreType$cuttag
    sepImage <- scoreMapSeparatorImage(config)
    scores <- list()
    summaryScores <- if(isCutTagScore){
        NULL
    } else if(isSampleScore) {
        scores$quantile <- getSampleScores(atac_metadata, config, coord, scoreTypeName, "quantile")
        scores$quantile$stageType # yes, this is stageType_delta (was named this way in score_functions.R)
    } else {
        getGenomeScores(sourceId, scoreTypeName, binI, stgmQuantile = TRUE, hicQuantile = TRUE)
    }
    imager::imappend(c(

        # a single row representing the fully aggregated score results for a scoreType
        # this is the only row for genome gc and txn scoreTypes
        # this row is omitted for cuttag scoreTypes
        list(
            scoreMapHeaderImage(scoreType, config),
            sepImage
        ),
        if(!isCutTagScore) list(
            scoreMapRowImage(
                scoreTypeName, isCutTagScore, !isSampleScore && scoreType$log10,
                binI, b, 
                summaryScores, 
                getSeriesSummaryColors, config, "summary", 
                labelLabel  = if(isSampleScore) "delta" else scoreType$trackSummaryLabel,
                legendLabel = if(isSampleScore) "round - elong" else scoreType$trackLegendLabel
            ),
            sepImage
        ) else list(),
        if(isSampleScore && config$Aggregate_By != "none") {
            metadata <- if(isCutTagScore) cuttag_metadata else atac_metadata
            samplesFilter <- if(isCutTagScore) metadata$samples$antibody_target == scoreTypeName else TRUE
            seriesAggNames <- getSeriesAggNames(metadata, config, samplesFilter)
            isGcrz <- startsWith(scoreTypeName, "gcrz")
            if(isGcrz && config$GCRZ_As_Quantiles) { # give user the option to plot GCRZ scores as quantiles
                scores$data <- getSampleScores(metadata, config, coord, scoreTypeName, "quantile")
                colorFn <- getSeriesQuantileColors
                labelLabel <- "quantile"
            } else {
                scores$data <- getSampleScores(metadata, config, coord, scoreTypeName, "score", isCutTagScore)
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
                        scoreTypeName, isCutTagScore, scoreType$log10,
                        binI, b, 
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
                if(!isGcrz && !isCutTagScore) { # 
                    if(is.null(scores$quantile)) scores$quantile <- getSampleScores(metadata, config, coord, scoreTypeName, "quantile")
                    lapply(
                        seriesAggNames,
                        function(seriesName) scoreMapRowImage(
                            scoreTypeName, isCutTagScore, FALSE,
                            binI, b,  
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
