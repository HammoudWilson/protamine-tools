#----------------------------------------------------------------------
# server components for the regionExpansion widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
regionExpansionServer <- function(id, sourceId, dinuc) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'regionExpansion'
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
# selected region profile plot
#----------------------------------------------------------------------
selectedRegion <- reactiveVal(NULL)
setSelectedRegion <- function(coord, d) {
    dists <- rowSums((d$xy - matrix(c(coord$x, coord$y), nrow(d$xy), 2, byrow=TRUE))^2)
    rowI <- which.min(dists)
    selectedRegion(list(
        region = d$regions[rowI],
        color  = d$color[rowI]
    ))
}
regionProfilePlotBox <- staticPlotBoxServer(
    "regionProfilePlotBox",
    title = TRUE,
    create = function() {
        sourceId <- sourceId()
        req(sourceId)
        ai <- paTss_ab_initio(sourceId)
        d <- selectedRegion()
        req(ai, d)
        startSpinner(session, message = "rendering profile")
        nStages <- length(ai$stages)
        stage_cols <- if(dinuc$input$rpkm_scaling == "scaled") {
            paste(ai$stages, "scaled", sep = '_')
        } else {
            paste(ai$stages, "rpkm", sep = '_')
        }
        y <- unlist(d$region[, .SD, .SDcols = stage_cols])
        ylim <- if(dinuc$input$rpkm_scaling == "scaled") range(y) else c(0, max(y))
        ylim <- ylim * 1.05
        par(mar = titledMar)
        regionProfilePlotBox$initializeFrame(
            xlim = c(1, nStages),
            ylim = ylim,
            xlab = "",
            ylab = if(dinuc$input$rpkm_scaling == "scaled") "log2(stage RPKM / mean RPKM)" else "RPKM",
            xaxt = "n",
            title = paste0(d$region$chrom, ":", d$region$start0, "-", d$region$end1)
        )
        addStageXAxis(ai$stages, ylim)
        if(dinuc$input$rpkm_scaling == "scaled") abline(h = 0, col = CONSTANTS$plotlyColors$grey)
        abline(v = d$region$stage_mean, col = d$color, lty = 1)
        regionProfilePlotBox$addLines(
            x = 1:nStages,
            y = y,
            col = d$color,
            lwd = 2
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# selected region insert footprint plot
#----------------------------------------------------------------------
createRegionPlot <- function(settings, plotBox) {
    region <- selectedRegion()$region
    req(region)
    startSpinner(session, message = "rendering region footprint")
    padding_bp <- 250
    coord <- list(chromosome = region$chrom, start = region$start0 + 1 - padding_bp, end = region$end1 + padding_bp)
    coord$range <- c(coord$start, coord$end)
    coord$width <- coord$end - coord$start + 1
    app$browser$createBrowserPlot(
        regionI = 1, 
        pngFile = plotBox$pngFile,
        externalCoord = coord
    )$layout
}
regionPlotBox <- mdiInteractivePlotBoxServer(
    "regionPlotBox",
    defaults = list(
        Plot_Frame = list(
            Width_Inches = 8,
            Height_Inches = 4
        )
    ),
    create = function(...) createRegionPlot(..., regionPlotBox) # a function or reactive that creates the plot as a png file using settings and helpers
)

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
list(
    setSelectedRegion = setSelectedRegion,
    regionProfilePlotBox = regionProfilePlotBox,
    regionPlotBox = regionPlotBox
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
