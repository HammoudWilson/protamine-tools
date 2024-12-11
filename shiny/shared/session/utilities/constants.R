# limit values
gcLimits <- c(0.25, 0.65) # optimal for mm39 genome, only used by normalizeGC step (other inherit from collate/scores)

# data types
refTypes <- c('genome', 'spike_in')
insertTypes <- c('subnucleosomal', 'nucleosomal')
getTypedStages <- function(sourceId) unlist(paScores(sourceId)$stageTypes)
getStageTypesByStage <- function(sourceId, stages) unlist(paScores(sourceId)$reverseStageTypes[stages])

# score types metadata
scoreTypes <- list(
    genome = list(
        gc = list(
            name = "Bin Fraction GC",
            gbBiasDependent = FALSE,
            distUnit = 0.01,
            class = "baseComposition",
            lim = gcLimits
        ),
        rpkm = list(
            name = "Transcription RPKM",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "transcription",
            lim = c(0, 10)
        )
    ),
    sample = list(
        cpm = list(
            name = "Counts Per Million",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "coverage",
            lim = c(0, 2)
        ),
        gcrz = list(
            name = "GC-Residual Z-Score",
            gbBiasDependent = TRUE, # thus, cannot be assessed until GC bias is established in app
            distUnit = 0.1,
            class = "coverage",
            lim = c(-2, 2)
        ),
        snif = list(
            name = "Subnucleosomal Insert Fraction",
            gbBiasDependent = FALSE,
            distUnit = 0.01,
            class = "insertSize",
            lim = c(0, 1)
        ),
        nrll = list(
            name = "Subnucleosomal vs. Nucleosomal NRLL",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "insertSize",
            lim = c(-2, 2)
        )
    )
)
getScoreLevel <- function(scoreType){
         if(scoreType %in% names(scoreTypes$genome)) 'genome'
    else if(scoreType %in% names(scoreTypes$sample)) 'sample'
    else "NA"
}
getScoreType <- function(sourceId, scoreType){
    scoreLevel <- getScoreLevel(scoreType)
    paScores(sourceId)$scoreTypes[[scoreLevel]][[scoreType]]
}

# standardized color palettes
paColors <- list(
    # BROWN  = rgb(0.2, 0,   0),   # extreme gain    
    # RED    = rgb(0.9, 0,   0),   # 3, CN3, full gain
    # YELLOW = rgb(0.8, 0.8, 0),   # 2.5, mosaic state
    # GREEN  = rgb(0,   0.8, 0),   # 2, CN2, no CNV
    # CYAN   = rgb(0,   0.6, 1),   # 1.5, mosaic state
    # BLUE   = rgb(0,   0,   1),   # 1, CN1, full loss
    # PURPLE = rgb(0.4, 0,   0.4), # extreme loss
    # BLACK  = rgb(0,   0,   0)    # 0, nothing    
    RED     = rgb(0.9, 0,   0),
    GREY    = rgb(0.75, 0.75, 0.75),
    BLUE    = rgb(0,   0,   1)
)
stageColors <- c(
    CONSTANTS$plotlyColors$black,
    CONSTANTS$plotlyColors$blue,
    CONSTANTS$plotlyColors$orange,
    CONSTANTS$plotlyColors$green,
    CONSTANTS$plotlyColors$purple
)
stageTypeColors <- c(
    CONSTANTS$plotlyColors$red,
    CONSTANTS$plotlyColors$blue
)
getSampleColorsByStage <- function(allSamples, samples){
    stages <- allSamples[, unique(stage)]
    colors <- sapply(allSamples$stage, function(x) stageColors[which(stages == x)])
    names(colors) <- allSamples$sample_name
    colors[samples$sample_name]
}
getStageColors <- function(allSamples, samples){
    stages <- allSamples[, unique(stage)]
    colors <- stageColors[1:length(stages)]
    names(colors) <- stages
    colors[samples[, unique(stage)]]
}
getStageTypeColors <- function(sourceId, allSamples, samples){
    allStageTypes <- getStageTypesByStage(sourceId, allSamples[, unique(stage)])
    stageTypes    <- getStageTypesByStage(sourceId,    samples[, unique(stage)])
    colors <- stageTypeColors[1:length(allStageTypes)]
    names(colors) <- allStageTypes
    colors[stageTypes]
}
nTrackMapColorsPerSide <- 30
trackMapColors <- list(
    low  = colorRampPalette(c(paColors$GREY, paColors$BLUE))(nTrackMapColorsPerSide + 1), # blue color is cold/depleted,
    high = colorRampPalette(c(paColors$GREY, paColors$RED))( nTrackMapColorsPerSide + 1) # red  color is hot/ enriched
)
