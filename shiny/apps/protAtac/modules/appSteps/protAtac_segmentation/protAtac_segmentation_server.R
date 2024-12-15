#----------------------------------------------------------------------
# server components for the protAtac_segmentation appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
protAtac_segmentationServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'protAtac_segmentation'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
trainingRegions <- reactiveValues()
sourceTrainingRegions <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    trainingRegions[[sourceId]]
})
observeEvent(input$clearTrainingRegions, {
    showUserDialog(
        "Clear All Training Regions?", 
        tags$p("Click OK to clear the training regions list."), 
        callback = function(...) {
            sourceId <- sourceId()
            req(sourceId)
            trainingRegions[[sourceId]] <- NULL
        },
        fade = FALSE
    )
})

#----------------------------------------------------------------------
# samples table
#----------------------------------------------------------------------
trainingRegionsTable <- bufferedTableServer(
    "trainingRegions",
    id,
    input,
    tableData = sourceTrainingRegions,
    selection = 'none',
    options = list(
        searching = FALSE  
    )
)

#----------------------------------------------------------------------
# appStep outcomes, saved to disk
#----------------------------------------------------------------------
trainingFileName <- "regionTraining.rds"
getTrainingFilePath <- function(sourceId) expandSourceFilePath(sourceId, trainingFileName)
trainingFilePath <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    getTrainingFilePath(sourceId)
})
invalidateTrainingResults <- reactiveVal(1)
getTrainingResults <- function(trainingFilePath = NULL, sourceId = NULL){
    if(is.null(trainingFilePath)) trainingFilePath <- getTrainingFilePath(sourceId)
    if(file.exists(trainingFilePath)) readRDS(trainingFilePath) else list()
}
writeTrainingResults <- function(sourceId, trainingResults){
    saveRDS(trainingResults, file = getTrainingFilePath(sourceId))
}
trainingResults <- reactive({ # all trainging results of all types
    invalidateTrainingResults()
    getTrainingResults( trainingFilePath = trainingFilePath() )
})
# trainingData <- reactive({ # the data for one type of training
#     trainingResults <- trainingResults()
#     sample <- sample()
#     trainingResults[[sample$sample_name]]
# })

#----------------------------------------------------------------------
# parse the state training distribtuions
#----------------------------------------------------------------------
getGenomeDistribution <- function(sourceId, scoreTypeName){
    getStageTypeDeltaScores(sourceId, scoreTypeName)[[1]]$dist
}
getTrainingDistribution <- function(sourceId, scoreTypeName, samples, isTrainingBin){
    stageTypeScores <- getStageTypeScores(sourceId, scoreTypeName, samples)
    delta <- stageTypeScores$round$score[isTrainingBin] - 
             stageTypeScores$elong$score[isTrainingBin]
    distUnit <- getScoreType(sourceId, scoreTypeName)$distUnit
    x <- data.table(x = floor(delta / distUnit) * distUnit)[, .(y = .N), keyby = .(x)]
    x[, y := y / sum(y)]
    x
}
trainingDistributions <- reactive({
    trainingRegions <- sourceTrainingRegions()
    req(trainingRegions, nrow(trainingRegions) > 0)
    sourceId <- sourceId()
    req(sourceId)
    bd <- paBinData(sourceId)
    setkeyv(trainingRegions, c("chrom", "start0", "end1"))
    binOverlaps <- foverlaps(
        bd$bins$genome,
        trainingRegions,
        type = "any",
        which = TRUE,
        mult = "first"
    )
    isTrainingBin <- !is.na(binOverlaps)
    list(
        gcrz = list(
            genome   = getGenomeDistribution(sourceId,   "gcrz"),
            training = getTrainingDistribution(sourceId, "gcrz", bd$samples, isTrainingBin)
        ),
        iisf = list(
            genome   = getGenomeDistribution(sourceId,   "iisf"),
            training = getTrainingDistribution(sourceId, "iisf", bd$samples, isTrainingBin)
        ),
        nrll = list(
            genome   = getGenomeDistribution(sourceId,   "nrll"),
            training = getTrainingDistribution(sourceId, "nrll", bd$samples, isTrainingBin)
        )
    )
})

#----------------------------------------------------------------------
# training region pre-HMM plots
#----------------------------------------------------------------------
distributionsPlot <- function(scoreTypeName, xlim) {
    plot <- staticPlotBoxServer(
        paste0(scoreTypeName, "Distributions"),
        maxHeight = "400px",
        points  = TRUE,
        margins = TRUE,
        title   = TRUE,
        create = function(){
            sourceId <- sourceId()
            req(sourceId)
            scoreType <- getScoreType(sourceId, scoreTypeName)
            d <- trainingDistributions()
            maxY <- max(
                d[[scoreTypeName]]$genome$y,
                d[[scoreTypeName]]$training$y,
                na.rm = TRUE
            )
            plot$initializeFrame(
                xlim = xlim,
                ylim = c(0, maxY * 1.05),
                xlab = paste(scoreType$label, scoreType$unit, "Delta"), 
                ylab = "Frequency"
            )
            plot$addLines( # addLines follows the same pattern, etc.
                x = d[[scoreTypeName]]$genome,
                pch = 16,
                cex = 0.2,
                col = paColors$BLUE
            )
            plot$addLines(
                x = d[[scoreTypeName]]$training,
                pch = 16,
                cex = 0.2,
                col = paColors$RED
            )
        }
    )
    plot
}
gcrzDistributionsPlot <- distributionsPlot("gcrz", c(-3, 5))
iisfDistributionsPlot <- distributionsPlot("iisf", c(-0.5, 0.25)) 
nrllDistributionsPlot <- distributionsPlot("nrll", c(-1.5, 0.5))

#----------------------------------------------------------------------
# use the distributions to segment the genome for two states that do and do not match the training regions
#----------------------------------------------------------------------
getLogEmissProbs <- function(sourceId, scoreTypeName, modelType){
    dists <- trainingDistributions()
    dist <- dists[[scoreTypeName]][[modelType]]
    minProb <- 1e-3
    dist[, y := pmax(minProb, y)] # ensure a minimum non-zero probability for all delta values for all models
    dist[, y := y / sum(y)]
    epdf <- approxfun(
        dist,
        yleft  = minProb, 
        yright = minProb
    )
    delta <- getStageTypeDeltaScores(sourceId, scoreTypeName)[[1]]$score
    log10(epdf(delta))
}
observeEvent(input$solveHMM, {
    req(trimws(input$trainingSetName))
    showUserDialog(
        "Solve HMM on Training Regions?", 
        tags$p("Click OK to segment the genome based on the training regions list."), 
        tags$p(paste("Training set: ", input$trainingSetName)),
        tags$p("The job will run ansynchronously for several minutes."),
        callback = function(...) {
            removeModal()
            solveTrainingHMM()
        },
        fade = FALSE
    )
})
solveTrainingHMM_monitor <- reactiveVal(NULL)
solveTrainingHMM_async <- function(sourceId, trainingSetName, hmm){
    list(
        sourceId = sourceId, 
        trainingSetName = trainingSetName,
        states = keyedViterbi(hmm) # 0 (genome) or 1(training)
    )
}
solveTrainingHMM <- function(){ # this is blocking at present, how long will it take?
    sourceId <- sourceId()
    req(sourceId)
    bd <- paBinData(sourceId)
    d <- trainingDistributions()
    hmm <- new_hmmEPTable(
        emissProbs = sapply(c("genome","training"), function(modelType){
            getLogEmissProbs(sourceId, "gcrz", modelType) + 
            getLogEmissProbs(sourceId, "iisf", modelType) +
            getLogEmissProbs(sourceId, "nrll", modelType)
        }), 
        transProb = as.double(settings$get("Segmentation","Transition_Probability")), 
        keys = bd$bins$genome[, chrom]
    )
    mdi_async(
        solveTrainingHMM_async,
        solveTrainingHMM_monitor,
        name = "solveTrainingHMM",
        header = TRUE,
        autoClear = NULL, # if not NULL, successful (but not error) header status icon is cleared after autoClear milliseconds
        sourceId = sourceId,
        trainingSetName = trimws(input$trainingSetName),
        hmm = hmm
    )
}
observeEvent(solveTrainingHMM_monitor(), {
    new <- solveTrainingHMM_monitor()
    req(!new$pending)
    if(new$success){
        new <- new$value

        dmsg("solveTrainingHMM_monitor")
        dstr(new)

        tr <- getTrainingResults(sourceId = new$sourceId)
        tr[[new$trainingSetName]] <- list(
            states = new$states,
            segments = paBinData(sourceId)$bins$genome[, .(
                chrom, 
                start0, 
                end1, 
                binI, 
                state = new$states
            )][, .(
                start0 = start0[1], 
                end1   = end1[.N], 
                state  = state[1]
            ), by = .(chrom, state)]
        )
        writeTrainingResults(new$sourceId, tr)
        invalidateTrainingResults(invalidateTrainingResults() + 1)
    } else {
        showUserDialog(
            "HMM Training Error", 
            tags$p("The HMM training job reported the following error:"),
            tags$p(new$message),
            type = "okOnly",
            size = if(nchar(new$message) > 200) "m" else "s"
        )
    }
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
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
    outcomes = list(),
    addTrainingRegion = function(sourceId, chrom, start0, end1, size_kb){
        tr <- trainingRegions[[sourceId]]
        if(is.null(tr)) tr <- data.table(
            chrom   = character(),
            start0  = integer(),
            end1    = integer(),
            size_kb = numeric()
        )
        trainingRegions[[sourceId]] <- rbind(
            tr,
            data.table(
                chrom   = chrom,
                start0  = start0,
                end1    = end1,
                size_kb = size_kb
            )
        )
    },
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
