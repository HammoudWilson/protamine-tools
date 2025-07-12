#----------------------------------------------------------------------
# server components for the intervalExpansion widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
intervalExpansionServer <- function(id, peaks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'intervalExpansion'
# settings <- activateMdiHeaderLinks( # uncomment as needed
#     session,
#     url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
#     dir = getAppStepDir(module), # for terminal emulator
#     envir = environment(), # for R console
#     baseDirs = getAppStepDir(module), # for code viewer/editor
#     settings = id, # for step-level settings
#     immediate = TRUE # plus any other arguments passed to settingsServer()
# )

#----------------------------------------------------------------------
# selected interval profile plot
#----------------------------------------------------------------------
selectedInterval <- reactiveVal(NULL)
setSelectedInterval <- function(coord, d) {
    dists <- rowSums((d$xy - matrix(c(coord$x, coord$y), nrow(d$xy), 2, byrow=TRUE))^2)
    rowI <- which.min(dists)
    selectedInterval(list(
        interval = d$intervals[rowI],
        color    = d$color[rowI]
    ))
}
profilePlotBox <- staticPlotBoxServer(
    "profilePlotBox",
    title = TRUE,
    create = function() {
        d <- selectedInterval()
        req(d)
        startSpinner(session, message = "rendering profile")
        stages <- peaks$stages()
        nStages <- length(stages)
        stage_cols <- if(peaks$input$rpkm_scaling == "scaled") {
            paste(stages, "scaled", sep = '_')
        } else {
            paste(stages, "rpkm", sep = '_')
        }
        y <- unlist(d$interval[, .SD, .SDcols = stage_cols])
        ylim <- if(peaks$input$rpkm_scaling == "scaled") range(y) else c(0, max(y))
        ylim <- ylim * 1.05
        par(mar = titledMar)
        profilePlotBox$initializeFrame(
            xlim = c(1, nStages),
            ylim = ylim,
            xlab = "",
            ylab = if(peaks$input$rpkm_scaling == "scaled") "log2(stage RPKM / mean RPKM)" else "RPKM",
            xaxt = "n",
            title = paste0(d$interval$chrom, ":", d$interval$start0, "-", d$interval$end1)
        )
        addStageXAxis(stages, ylim)
        if(peaks$input$rpkm_scaling == "scaled") abline(h = 0, col = CONSTANTS$plotlyColors$grey)
        abline(v = d$interval$stage_mean, col = d$color, lty = 1)
        profilePlotBox$addLines(
            x = 1:nStages,
            y = y,
            col = d$color,
            lwd = 2
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# selected interval insert footprint plot
#----------------------------------------------------------------------
browserPlotSettings <- list(
    Expansion = list(
        Padding_bp = list(
            type = "numericInput",
            value = 1000,
            min = 0,
            max = 10000,
            step = 100
        )
    )
)
createBrowserPlot <- function(settings, plotBox) {
    interval <- selectedInterval()$interval
    req(interval)
    startSpinner(session, message = "rendering embedded browser")
    padding_bp <- browserPlotBox$settings$get("Expansion", "Padding_bp")
    coord <- list(chromosome = interval$chrom, start = interval$start0 + 1 - padding_bp, end = interval$end1 + padding_bp)
    coord$range <- c(coord$start, coord$end)
    coord$width <- coord$end - coord$start + 1
    app$browser$createBrowserPlot(
        regionI = 1, 
        pngFile = plotBox$pngFile,
        externalCoord = coord
    )$layout
}
browserPlotBox <- mdiInteractivePlotBoxServer(
    "browserPlotBox",
    settings = browserPlotSettings,
    defaults = list(
        Plot_Frame = list(
            Width_Inches = 10,
            Height_Inches = 4
        )
    ),
    create = function(...) createBrowserPlot(..., browserPlotBox) # a function or reactive that creates the plot as a png file using settings and helpers
)

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
list(
    setSelectedInterval = setSelectedInterval,
    profilePlotBox = profilePlotBox,
    browserPlotBox = browserPlotBox
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
