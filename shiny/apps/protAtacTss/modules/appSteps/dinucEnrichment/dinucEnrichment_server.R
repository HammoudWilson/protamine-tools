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
dinuc <- dinucRegionsSelectorBoxServer("dinucRegions", sourceId)
regionsBedTable <- regionsBedTableServer("regionsBedTable", sourceId)
clusterProfilePlotBox <- clusterProfilePlotBoxServer("clusterProfilePlotBox", sourceId, dinuc)
regionExpansion <- regionExpansionServer("regionExpansion", sourceId, dinuc)

#----------------------------------------------------------------------
# overlap (fraction) enrichment plot
#----------------------------------------------------------------------
overlapPlotSettings <- list(
    Enrichment = list(
        Enrichment_X_Axis = list(
            type = "selectInput",
            choices = c("quantile", "stage_mean"),
            value = "stage_mean"
        )
    )
)
overlapPlotData <- reactive({
    d <- clusterProfilePlotBox$data()
    bedData <- regionsBedTable$data()
    req(d, bedData)
    startSpinner(session, message = "getting overlaps")
    overlaps <- foverlaps(
        d$regions,        # larger table, smaller intervals
        bedData,          # keyed, smaller table, larger intervals
        type = "any",     # detect any overlap
        mult = "first",   # we only need to know if there is at least one overlap
        nomatch = NA,     # keep non-matching rows
        which = TRUE      # return indices instead of joined data
    )
    d$regions[, hasOverlap := !is.na(overlaps)] 
    d$fractions <- d$regions[, .(
        frac = sum(hasOverlap) / .N,
        stage_mean = mean(stage_mean, na.rm = TRUE)
    ), keyby = c(d$clusterCol)]
    d
})
overlapPlotBox <- staticPlotBoxServer(
    "overlapPlotBox",
    settings = c(overlapPlotSettings),
    create = function() {
        d <- overlapPlotData()
        bedMetadata <- regionsBedTable$metadata()
        req(d, bedMetadata)
        startSpinner(session, message = "rendering overlap")
        isQuantile <- overlapPlotBox$settings$get("Enrichment", "Enrichment_X_Axis") == "quantile"
        par(mar = titledMar)
        overlapPlotBox$initializeFrame(
            xlim = if(isQuantile) c(0.5, d$nClusters + 0.5) else c(1, d$nStages),
            ylim = c(0, max(d$fractions$frac, na.rm = TRUE) * 1.05),
            xlab = if(isQuantile) d$clusterType else "Stage Mean",
            ylab = "Fraction Overlapping BED"
        )
        abline(h = bedMetadata$fraction_of_genome, col = CONSTANTS$plotlyColors$grey, lty = 2)
        x <- if(isQuantile) 1:d$nClusters else d$fractions$stage_mean
        y <- d$fractions$frac
        overlapPlotBox$addLines(
            x = x,
            y = y,
            col = CONSTANTS$plotlyColors$grey
        )
        overlapPlotBox$addPoints(
            x = x,
            y = y,
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
        Missing_Overlap_Score = list(
            type = "selectInput",
            choices = c("zero", "omit"),
            value = "omit"
        ),
        Invert_Overlaps = list(
            type = "checkboxInput",
            value = FALSE
        )
    )
)
scorePlotData <- reactive({
    d <- clusterProfilePlotBox$data()
    bedData <- regionsBedTable$data()
    req(d, bedData, ncol(bedData) >= 5)
    startSpinner(session, message = "getting scores")
    regionCols <- c("chrom", "start0", "end1", "stage_mean", d$clusterCol)
    d$scoreCol <- names(bedData)[5]
    overlaps <- foverlaps(
        d$regions[, ..regionCols], # larger table, smaller intervals
        bedData[, c(1:3, 5)],      # keyed, smaller table, larger intervals
        type = "any",  # detect any overlap
        mult = "all",  # multiple span matches aggregated to a single value below
        nomatch = NA,  # keep non-matching rows (for now)
        which = FALSE  # return a joined data frame
    )
    naScore <- if(scorePlotBox$settings$get("Enrichment", "Missing_Overlap_Score") == "zero") 0.0 else NA_real_
    d$scores <- overlaps[, {
        hasOverlap <- !is.na(start0[1]) && (
            dinuc$input$allow_partial_overlap ||
            any(i.start0 >= start0 & i.end1 <= end1)
        )
        .(
            score = if(hasOverlap) {
                ss <- if(dinuc$input$allow_partial_overlap) .SD[[d$scoreCol]] else mapply(function(s0, e1, s) {
                    if(i.start0 >= s0 & i.end1 <= e1) s else NA_real_
                }, start0, end1, .SD[[d$scoreCol]])
                max(ss, na.rm = TRUE) 
            } else naScore,
            hasOverlap = hasOverlap,
            cluster    = .SD[[d$clusterCol]][1],
            stage_mean = stage_mean[1]
        )
    }, keyby = .(chrom, i.start0, i.end1)]
    if(scorePlotBox$settings$get("Enrichment", "Invert_Overlaps")){
        d$scores[, score := ifelse(hasOverlap, NA_real_, runif(.N, 0, 1))]
    }
    d$regions <- d$regions[order(stage_mean)]
    d$scores  <- d$scores[ order(stage_mean)]
    d$scores[, color := paRainbow(.N, 0.5)]
    d$medians <- d$scores[, .(
        stage_mean = median(stage_mean, na.rm = TRUE),
        score      = median(score, na.rm = TRUE)
    ), keyby = cluster]
    d$xy <- as.matrix(d$scores[, .(stage_mean, score)])
    d$color <- d$scores$color
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
            ylim = range(d$scores$score, na.rm = TRUE) * 1.05,
            xlab = "Stage Mean",
            ylab = gsub("_", " ", d$scoreCol),
            title = paste(
                format(sum(d$scores$hasOverlap), big.mark = ","), 
                "of", 
                format(nrow(d$scores), big.mark = ","),
                "dinucs"
            )
        )
        scorePlotBox$addPoints(
            x   = d$scores$stage_mean,
            y   = d$scores$score,
            col = d$scores$color
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
    regionExpansion$setSelectedRegion(scorePlotBox$plot$click()$coord, scorePlotData())
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
        regionExpansion$regionProfilePlotBox$settings$replace(bm$outcomes$regionProfilePlotBoxSettings)
        regionExpansion$regionPlotBox$settings$replace(bm$outcomes$regionPlotBoxSettings)
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
        regionProfilePlotBoxSettings = regionExpansion$regionProfilePlotBox$settings$all_,
        regionPlotBoxSettings = regionExpansion$regionPlotBox$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
