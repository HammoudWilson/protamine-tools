# limit values
gcLimits <- c(0.25, 0.65) # optimal for mm39 genome, only used by normalizeGC step (other inherit from collate/scores)

# data types
refTypes <- c('genome', 'spike_in')
insertTypes <- c('subnucleosomal', 'nucleosomal')

# score types metadata
scoreTypes <- list(
    genome = list(
        gc = list(
            name = "Bin Fraction GC",
            gbBiasDependent = FALSE,
            distUnit = 0.005,
            class = "baseComposition",
            valueLim = gcLimits,
            normValue = "z"
        ),
        rpkm = list(
            name = "Transcription RPKM",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "transcription",
            valueLim = c(0, 10),
            normValue = "quantile"
        )
    ),
    sample = list(
        cpm = list(
            name = "Counts Per Million",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "coverage",
            valueLim = c(0, 1.5),
            deltaLim = c(-0.5, 0.5),
            normValue = "quantile"
        ),
        gcrz = list(
            name = "GC-Residual Z-Score",
            gbBiasDependent = TRUE, # thus, cannot be assessed until GC bias is established in app
            distUnit = 0.1,
            class = "coverage",
            valueLim = c(-3, 3),
            deltaLim = c(-3, 3),
            normValue = "z"
        ),
        snif = list(
            name = "Subnucleosomal Insert Fraction",
            gbBiasDependent = FALSE,
            distUnit = 0.01,
            class = "insertSize",
            valueLim = c(0, 1),
            deltaLim = c(-1, 0.5),
            normValue = "z"
        ),
        nrll = list(
            name = "Subnucleosomal vs. Nucleosomal NRLL",
            gbBiasDependent = FALSE,
            distUnit = 0.1,
            class = "insertSize",
            valueLim = c(-1.5, 1.5),
            deltaLim = c(-1.5, 0),
            normValue = "z"
        )
    )
)


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
nTrackMapColorsPerSide <- 30
trackMapColors <- list(
    low  = colorRampPalette(c(paColors$GREY, paColors$BLUE))(nTrackMapColorsPerSide + 1), # blue color is cold/depleted,
    high = colorRampPalette(c(paColors$GREY, paColors$RED))( nTrackMapColorsPerSide + 1) # red  color is hot/ enriched
)
