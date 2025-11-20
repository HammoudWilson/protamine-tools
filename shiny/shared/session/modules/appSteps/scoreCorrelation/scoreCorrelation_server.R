#----------------------------------------------------------------------
# server components for the scoreCorrelation appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
scoreCorrelationServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'scoreCorrelation'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    size = "m"
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data sources
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "dataSourceTable", 
    selection = "single"
)
nBinsPerPoint <- reactive({
    as.integer(settings$get("Score_Correlation", "N_Bins_Per_Point"))
})
nPoints <- reactive({
    settings$get("Score_Correlation", "N_Plot_Points")
})
minQuantile <- reactive({
    settings$get("Score_Correlation", "Data_Range_Min_Quantile")
})
randomBinI <- reactive({
    sourceId <- req(sourceId())
    req(sourceId)
    bins <- paScores_bins(sourceId)
    includedAutosomeBins <- getIncludedAutosomeBins_scores(bins)
    bins[, isIncluded := includedAutosomeBins & !is.na(stgm)]
    nBinsPerPoint <- nBinsPerPoint()
    nPoints <- nPoints()
    binI <- if(nBinsPerPoint == 1){
        bins[, which(isIncluded)]
    } else {
        bins[, 
            .(
                binI = if(.N < nBinsPerPoint) NA_real_
                       else seq(.I[1] + nBinsPerPoint - 1, .I[.N], nBinsPerPoint) - nBinsPerPoint + 1
            ), 
            by = .(chrom, rleid(isIncluded))
        ][!is.na(binI), as.integer(binI)]
    }
    sample(binI, min(nPoints, length(binI)))
})
binIndices <- reactive({
    nBinsPerPoint <- nBinsPerPoint()
    randomBinI <- req(randomBinI())
    outer(randomBinI, 0:(nBinsPerPoint - 1), `+`)
})

#----------------------------------------------------------------------
# plot outputs
#----------------------------------------------------------------------
roundMinusElong <- "round - elong"
getScoreTypeData <- function(scoreTypeName, stage_, isZAxis = FALSE){
    if(isZAxis && scoreTypeName == "gc") scoreTypeName <- "gc_z"
    sourceId <- req(sourceId())
    nBinsPerPoint <- nBinsPerPoint()
    randomBinI <- req(randomBinI())
    scoreLevel <- getScoreLevel(scoreTypeName)
    scoreType <- scoreTypes[[scoreLevel]][[scoreTypeName]]
    isSampleScore <- scoreLevel == "sample"
    isCutTagScore <- !is.null(scoreType$cuttag) && scoreType$cuttag
    isDeltaScore <- isSampleScore && stage_ == roundMinusElong
    scores <- if(isSampleScore) {
        if(isDeltaScore){
            req(!isCutTagScore)
            getStageTypeDelta_allBins(sourceId, scoreTypeName)
        } else {
            getSampleScores_allBins(sourceId, scoreTypeName, stage_)
        }
    } else {
        # do not apply quantile transform to genome scores for X and Y axes, but do for Z axis coloring
        getGenomeScores(sourceId, scoreTypeName, TRUE, stgmQuantile = isZAxis, hicQuantile = isZAxis) 
    }
    if(isCutTagScore){
        atac_metadata <- paScores_metadata(sourceId)
        samples <- paCutTag_samples(sourceId)
        k_bin <- atac_metadata$env$BIN_SIZE / 1e3
        m_reads <- samples[antibody_target == scoreTypeName & stage == stage_, sum(n_reads) / 1e6]
    }
    scores <- if(nBinsPerPoint == 1){
        if(isCutTagScore) scores[randomBinI] / k_bin / m_reads else scores[randomBinI]
    } else {
        binIndices <- binIndices()
        nPoints <- nrow(binIndices) # may be less than nPoints() due to insufficient aggregated bins
        isLog10 <- !is.null(scoreType$log10) && scoreType$log10
        if(isCutTagScore) {
            rowMeans(matrix(scores[binIndices] / k_bin / m_reads, nrow = nPoints))
        } else if(isLog10) log10(rowMeans(matrix(10**scores[binIndices], nrow = nPoints)))
          else                   rowMeans(matrix(    scores[binIndices], nrow = nPoints))
    }
    list(
        scoreTypeName = scoreTypeName, 
        scoreType     = getScoreType(scoreTypeName),
        scoreLevel    = scoreLevel,
        isSampleScore = isSampleScore,
        isCutTagScore = isCutTagScore,
        isDeltaScore  = isDeltaScore,
        scores        = scores
    )
}
getScoreTypeAxis <- function(scoreTypeName, stage_, isZAxis = FALSE){
    d <- getScoreTypeData(scoreTypeName, stage_, isZAxis)
    unit <- d$scoreType[if(!is.null(d$scoreType$corrUnit)) "corrUnit" else "unit"]
    label <- paste(d$scoreType$label, unit)
    minQuantile <- minQuantile()
    min <- quantile(d$scores,     minQuantile, na.rm = TRUE)
    max <- quantile(d$scores, 1 - minQuantile, na.rm = TRUE)
    width <- max - min
    min <- min - 0.05 * width
    max <- max + 0.05 * width
    c(d, list(
        unit          = unit,
        label         = if(d$isSampleScore) {
            if(d$isDeltaScore) paste(label, "(Round - Elong)")
                          else paste(stage_, label)
        } else label,
        lim           = c(min, max),
        col = if(isZAxis) {
            config <- list(Max_Z_Score = 2, Min_Txn_Log10_CPM = -3, Max_Txn_Log10_CPM = 3)
            if(d$isSampleScore) getSeriesSampleColors(d$scoreTypeName, d$scores, config) else 
                               getSeriesSummaryColors(d$scoreTypeName, d$scores, config)
        } else NULL
    ))
}

#----------------------------------------------------------------------
# correlation plot output
#----------------------------------------------------------------------
xAxis <- reactive({
    getScoreTypeAxis(input$xScoreType, input$xScoreStage)
})
yAxis <- reactive({
    getScoreTypeAxis(input$yScoreType, input$yScoreStage)
})
zAxis <- reactive({
    getScoreTypeAxis(input$zScoreType, input$zScoreStage, TRUE)
})
getAxisLimits <- function(axis){
    limType <- settings$get("Score_Correlation", "Set_Axis_Limits")
    allowCorrRange <- is.null(axis$scoreType$allowCorrRange) || axis$scoreType$allowCorrRange
    if(!allowCorrRange || limType == "from_defaults"){
        corrLim <- if(!is.null(axis$scoreType$corrLim)) "corrLim" else "valueLim"
        axis$scoreType[[if(axis$isDeltaScore) "deltaLim" else corrLim]]
    } else {
        axis$lim
    }
}
zAxisTitle <- function(axis){
    if(is.null(axis$scoreType$trackLegendLabel)) axis$scoreType$trackHeaderLabel else paste(
        axis$scoreType$trackHeaderLabel, 
        axis$scoreType$trackLegendLabel, 
        sep = ", "
    )
}
plotCorrelation <- function(plot, xAxis, yAxis, zAxis){
    startSpinner(session, message = "plotting correlation")
    xlim <- getAxisLimits(xAxis)
    ylim <- getAxisLimits(yAxis)
    plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = xAxis$label,
        ylab = yAxis$label,
        title = zAxisTitle(zAxis)
    )
    plot$addPoints(
        x = xAxis$scores,
        y = yAxis$scores,
        # col = CONSTANTS$plotlyColors$blue %>% addAlphaToColor(0.1),
        col = zAxis$col,
        pch = 16,
        cex = 0.5
    )

    # overplot a unity line when comparing samples of the same score type
    if(
        (xAxis$isCutTagScore && yAxis$isCutTagScore) || 
         xAxis$scoreTypeName == yAxis$scoreTypeName
    ){
        plot$addLines(
            x = xAxis$lim,
            y = xAxis$lim,
            col = CONSTANTS$plotlyColors$grey,
            lwd = 1.5,
            lty = 1
        )
    }

    # overplot a linear model fit
    fit <- lm(y ~ x, data.frame(x = xAxis$scores, y = yAxis$scores))
    xFit <- seq(    
        from = min(xAxis$scores, na.rm = TRUE), 
        to   = max(xAxis$scores, na.rm = TRUE), 
        length.out = 100
    )
    yFit <- predict(fit, newdata = data.frame(x = xFit))
    plot$addLines(
        x = xFit,
        y = yFit,
        col = CONSTANTS$plotlyColors$black,
        lwd = 1
    )
    stopSpinner(session)
}
correlationPlot <- staticPlotBoxServer(
    "correlationPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        plotCorrelation(
            plot     = correlationPlot, 
            xAxis    = xAxis(),
            yAxis    = yAxis(),
            zAxis    = zAxis()
        )
    }
)

#----------------------------------------------------------------------
# PCA plot output
#----------------------------------------------------------------------
pca <- reactive({
    sourceId <- req(sourceId())
    req(input$pcaScoreTypes, length(input$pcaScoreTypes) >= 2)
    m <- do.call(cbind, lapply(input$pcaScoreTypes, function(scoreTypeName){
        startSpinner(session, message = paste("getting", scoreTypeName, "for PCA"))
        scoreLevel <- getScoreLevel(scoreTypeName)
        scoreType <- scoreTypes[[scoreLevel]][[scoreTypeName]]
        isCutTagScore <- !is.null(scoreType$cuttag) && scoreType$cuttag
        if(isCutTagScore){
            # getScoreTypeData(scoreTypeName, "earliest_ES")$scores -
            # getScoreTypeData(scoreTypeName, "early_ES")$scores
            cbind(
                getScoreTypeData(scoreTypeName, "earliest_ES")$scores,
                getScoreTypeData(scoreTypeName, "early_ES")$scores
            )
        } else if(scoreTypeName == "gcrz_obs") {
            # TODO: include all stages as features?
            getScoreTypeData(scoreTypeName, roundMinusElong)$scores
        } else {
            getScoreTypeData(scoreTypeName, "NA")$scores
        }
    }))
    nBefore <- nrow(m)
    keptRows <- apply(m, 1, function(x) all(is.finite(x)))
    m <- m[keptRows,]
    message(paste("removed", nBefore - nrow(m), "of", nBefore, "rows with non-finite values for PCA"))
    pca <- prcomp(m, retx = TRUE, center = TRUE, scale. = TRUE)
    pca$keptRows <- keptRows
    pca
})
plotPca <- function(plot, pca, zAxis){
    I <- sample(1:nrow(pca$x), nrow(pca$x))
    xlim <- range(pca$x[, 1], na.rm = TRUE)
    ylim <- range(pca$x[, 2], na.rm = TRUE)
    startSpinner(session, message = "plotting PCA")
    plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = "Principal Component 1",
        ylab = "Principal Component 2",
        title = zAxisTitle(zAxis)
    )
    plot$addPoints(
        x = pca$x[I, 1],
        y = pca$x[I, 2],
        col = zAxis$col[pca$keptRows][I],
        pch = 16,
        cex = 0.5
    )
    stopSpinner(session)
}
pcaPlot <- staticPlotBoxServer(
    "pcaPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = TRUE,
    title   = TRUE,
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        plotPca(
            plot  = pcaPlot,
            pca   = pca(),
            zAxis = zAxis()
        )
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        correlationPlot$settings$replace(bm$outcomes$correlationPlotSettings)
        pcaPlot$settings$replace(bm$outcomes$pcaPlotSettings)
    }
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(
        correlationPlotSettings = correlationPlot$settings$all_,
        pcaPlotSettings         = pcaPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
