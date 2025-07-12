#----------------------------------------------------------------------
# server components for the dinucEnrichment appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
dinucEnrichmentServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'dinucEnrichment'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
chains   <- dinucChainsSelectorBoxServer("dinucChains", sourceId)
regionsBedTable       <- regionsBedTableServer("regionsBedTable", sourceId)
clusterProfilePlotBox <- clusterProfilePlotBoxServer("clusterProfilePlotBox", chains)
intervalExpansion     <- intervalExpansionServer("intervalExpansion", chains)

#----------------------------------------------------------------------
# overlap assessment of x (dinuc chains) vs. y (BED regions)
#----------------------------------------------------------------------
overlapType <- reactive({
    # when allow_partial_overlap is FALSE, x must be completely contained within y
    # i.e., dinuc chains must be completely contained within a single BED region
    if(regionsBedTable$input$allow_partial_overlap) "any"
    else "within"
})
unfilteredOverlaps <- reactive({
    intervals <- chains$unfilteredIntervals()
    bedData <- regionsBedTable$data()
    req(intervals, bedData)
    paTss_interval_overlaps(chains, bedData, intervals, overlapType())

    # startSpinner(session, message = "finding overlaps")
    # bed3Cols   <- c("chrom", "start0", "end1")
    # clusterCol <- if(chains$input$include_unpassed) "quantile_unfiltered" else "quantile_filtered"
    # intervalCols <- c(bed3Cols, "stage_mean", clusterCol, "indexI")
    # scoreCol <- names(bedData)[5]
    # bedCols    <- c(bed3Cols, scoreCol)
    # scores <- foverlaps(
    #     x = intervals[, ..intervalCols],
    #     y = bedData[, ..bedCols],
    #     type = overlapType(), # any or within, see above
    #     mult = "all",         # multiple y matches aggregated to a single value below
    #     nomatch = NA,         # keep non-matching rows (for now)
    #     which = FALSE         # return left outer-join of x and y
    # )[, 
    #     {
    #         hasOverlap <- !is.na(start0[1])
    #         .(
    #             score      = if(hasOverlap) max(.SD[[scoreCol]], na.rm = TRUE) else NA_real_,
    #             hasOverlap = hasOverlap,
    #             cluster    = .SD[[clusterCol]][1],
    #             stage_mean = stage_mean[1]
    #         )
    #     }, 
    #     keyby = .(chrom, i.start0, i.end1)
    # ]
    # list(
    #     scoreCol = scoreCol,
    #     scores = scores
    # )
})
overlaps <- reactive({
    d <- c(
        clusterProfilePlotBox$data(),
        unfilteredOverlaps()
    )
    req(d)
    startSpinner(session, message = "filtering overlaps")
    d$scores <- d$scores[chains$intervals()[, indexI]]
    d$overlapSummary <- paste(
        format(sum(d$scores$hasOverlap), big.mark = ","), 
        "of", 
        format(nrow(d$scores), big.mark = ","),
        "peaks with BED overlap"
    )
    d
})

#----------------------------------------------------------------------
# overlap (fraction) enrichment plot
#----------------------------------------------------------------------
overlapPlotSettings <- list(
    Enrichment = list(
        Enrichment_X_Axis = list(
            type = "selectInput",
            choices = c("stage_mean", "quantile"),
            value = "stage_mean"
        )
    )
)
overlapPlotData <- reactive({
    d <- overlaps()
    clusters <- d$scores[, .(
        stage_mean   = mean(stage_mean, na.rm = TRUE),
        frac_overlap = sum(hasOverlap) / .N
    ), keyby = .(cluster)]
    if(overlapPlotBox$settings$get("Enrichment", "Enrichment_X_Axis") == "quantile") {
        d$xy   <- as.matrix(clusters[, .(cluster,    frac_overlap)])
        d$xlim <- c(0.5, d$nClusters + 0.5)
        d$xlab <- d$clusterType
    } else {
        d$xy   <- as.matrix(clusters[, .(stage_mean, frac_overlap)])
        d$xlim <- c(1, d$nStages)
        d$xlab <- "Stage Mean"
    }
    d
})
overlapPlotBox <- staticPlotBoxServer(
    "overlapPlotBox",
    settings = c(overlapPlotSettings),
    create = function() {
        d <- overlapPlotData()
        bedMetadata <- regionsBedTable$metadata()
        req(d, bedMetadata)
        startSpinner(session, message = "rendering fractions")
        par(mar = titledMar)
        overlapPlotBox$initializeFrame(
            xlim = d$xlim,
            ylim = c(0, max(d$xy[, 2], bedMetadata$fraction_of_genome, na.rm = TRUE) * 1.05),
            xlab = d$xlab,
            ylab = "Fraction Overlapping BED",
            title = d$overlapSummary
        )
        abline(h = bedMetadata$fraction_of_genome, col = CONSTANTS$plotlyColors$grey, lty = 2)
        overlapPlotBox$addLines(
            x = d$xy[, 1],
            y = d$xy[, 2],
            col = CONSTANTS$plotlyColors$grey
        )
        overlapPlotBox$addPoints(
            x = d$xy[, 1],
            y = d$xy[, 2],
            col = d$colors,
            pch = 19,
            cex = 1.5
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# overlap (fraction) enrichment plot
#----------------------------------------------------------------------
scorePlotSettings <- list(
    Enrichment = list(
        Invert_Overlaps = list(
            type = "checkboxInput",
            value = FALSE
        )
    )
)
scorePlotData <- reactive({
    d <- overlaps()
    d$scores    <- d$scores[order(stage_mean)]
    d$intervals <- d$intervals[order(stage_mean)]
    if(scorePlotBox$settings$get("Enrichment", "Invert_Overlaps")) {
        d$scores[, score := ifelse(hasOverlap, NA_real_, runif(.N, 0, 1))]
        d$scoreCol <- "Random Score"
    }
    d$xy <- as.matrix(d$scores[, .(stage_mean, score)])
    d$color <- paRainbow(nrow(d$scores), 0.5)
    d$medians <- d$scores[, .(
        stage_mean = median(stage_mean, na.rm = TRUE),
        score      = median(score, na.rm = TRUE)
    ), keyby = .(cluster)]
    d
})
scorePlotBox <- mdiInteractivePlotBoxServer(
    "scorePlotBox",
    click = TRUE,
    points = TRUE,
    lines = TRUE,
    settings = scorePlotSettings,
    defaults = c(paInteractiveFrame, list(
        Points_and_Lines = list(
            Point_Type = 19,
            Point_Size = 0.75
        )
    )),
    create = function(...) {
        d <- scorePlotData()
        bedMetadata <- regionsBedTable$metadata()
        req(d, bedMetadata)
        startSpinner(session, message = "rendering scores")
        layout <- scorePlotBox$initializePng(
            mar = titledMar
        ) %>% scorePlotBox$initializeFrame(
            xlim = c(1, d$nStages),
            ylim = range(d$xy[, 2], na.rm = TRUE) * 1.05,
            xlab = "Stage Mean",
            ylab = gsub("_", " ", d$scoreCol),
            title = d$overlapSummary
        )
        randomOrder <- sample(1:nrow(d$xy))
        scorePlotBox$addPoints(
            x   = d$xy[randomOrder, 1],
            y   = d$xy[randomOrder, 2],
            col = d$color[randomOrder]
        )
        scorePlotBox$addLines(
            x   = d$medians$stage_mean,
            y   = d$medians$score,
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
        scorePlotBox$finishPng(layout)
    }
)
observeEvent(scorePlotBox$plot$click(), {
    intervalExpansion$setSelectedInterval(scorePlotBox$plot$click()$coord, scorePlotData())
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    if(!is.null(bm$outcomes)){
        clusterProfilePlotBox$settings$replace(bm$outcomes$clusterProfilePlotBoxSettings)
        scorePlotBox$settings$replace(bm$outcomes$scorePlotBoxSettings)
        intervalExpansion$profilePlotBox$settings$replace(bm$outcomes$profilePlotBoxSettings)
        intervalExpansion$browserPlotBox$settings$replace(bm$outcomes$browserPlotBoxSettings)
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
    # settings = settings$all_,
    outcomes = list(
        clusterProfilePlotBoxSettings = clusterProfilePlotBox$settings$all_,
        scorePlotBoxSettings = scorePlotBox$settings$all_,
        profilePlotBoxSettings = intervalExpansion$profilePlotBox$settings$all_,
        browserPlotBoxSettings = intervalExpansion$browserPlotBox$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
