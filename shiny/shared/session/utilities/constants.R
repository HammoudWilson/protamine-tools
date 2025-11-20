# limit values
gcLimits <- c(0.25, 0.6) # optimal for mm39 genome; for normalization plots only, all other inherit from pipeline env

# empirically determine insert size values
meanNucleosomeFragSize   <- 205
meanDinucleosomeFragSize <- 375

# fixed genome values
genome_metadata <- list(
    mm39 = list(
        total_bp = 2723431143,
        included_bp = 2542229700
    )
)

# data types
refTypes <- c('genome', 'spike_in')
insertTypes <- c('all_inserts', 'intermediate')

# score types metadata
scoreTypes <- list(
    genome = list(
        # see definition of gc_z below, copied from gc
        gc = list(
            label = "GC",
            log10 = FALSE,
            unit = "Percent",
            trackHeaderLabel = "GC Percent",
            trackSummaryLabel = "Z Score",
            trackLegendLabel = "red/blue=high/low GC",
            class = "baseComposition",
            valueLim = gcLimits,
            corrLim = c(0.3, 0.575),
            allowCorrRange = FALSE
        ),
        txn = list(
            label = "Round Spermatid PRO-seq",
            log10 = TRUE,
            unit = "log10 RPKM",
            trackHeaderLabel = "Round Spermatid PRO-seq",
            trackSummaryLabel = "log10 RPKM",
            trackLegendLabel = "red=transcribed",
            class = "transcription",
            valueLim = log10(c(1e-3, 1e1)),
            corrUnit = "(log10 RPKM)",
            allowCorrRange = FALSE
        ),
        hic = list(
            label = "HiC Compartment Score",
            log10 = FALSE,
            unit = "",
            trackHeaderLabel = "HiC Compartment Score",
            trackSummaryLabel = "quantile",
            trackLegendLabel = "red/blue=A/B",
            class = "compartment",
            valueLim = c(-0.1, 0.1),
            allowCorrRange = FALSE
        ),
        stgm = list(
            label = "ATAC Stage Mean",
            log10 = FALSE,
            unit = "",
            trackHeaderLabel = "ATAC Stage Mean",
            trackSummaryLabel = "quantile",
            trackLegendLabel = "red/blue=early/late",
            class = "timecourse",
            valueLim = c(1.5, 5)
        )
    ),
    sample = list(
        gcrz_obs = list(
            distUnit = 0.1,
            include = c("quantile"),
            log10 = FALSE,
            label = "ATAC GC Residual",
            unit = "Z Score",
            enrichmentLabel = "ATAC GC Regression (Z)",
            trackHeaderLabel = "ATAC GC Regression",
            trackScoreLabel = "Z Score",
            trackLegendLabel = "red/blue=high/low Z-score",
            class = "coverage",
            valueLim = c(-3, 3),
            deltaLim = c(-3, 3)
        ),
        gcrz_wgt = list(
            distUnit = 0.1,
            include = c("quantile"),
            log10 = FALSE,
            label = "ATAC GC Residual",
            unit = "Z Score",
            enrichmentLabel = "ATAC GC Regression Tn5 Weighted (Z)",
            trackHeaderLabel = "ATAC GC Regression Tn5 Weighted",
            trackScoreLabel = "Z Score",
            class = "coverage",
            valueLim = c(-3, 3),
            deltaLim = c(-3, 3)
        ),
        nrll = list(
            distUnit = 0.1,
            log10 = FALSE,
            label = "ATAC Insert Size", # "Protamine Transition Enrichment",
            unit = "NRLL",
            enrichmentLabel = "Protamine Transition (NRLL)",
            trackHeaderLabel = "ATAC Small Insert Enrichment",
            trackScoreLabel = "NRLL",
            class = "insertSize",
            valueLim = c(-1.5, 1.5),
            deltaLim = c(-1.1, 0.1),
            corrUnit = "(NRLL)"
        ),
        H2B = list(
            distUnit = 1,
            include = c("quantile"),
            log10 = FALSE,
            label = "H2B Cut&Tag",
            unit = "# Reads",
            enrichmentLabel = "H2B (log10 RPKM)",
            trackHeaderLabel = "H2B Cut&Tag",
            trackScoreLabel = "RPKM",
            class = "chromatin",
            valueLim = c(0,50),
            deltaLim = c(-3, 3),
            corrLim  = c(0, 2),
            corrUnit = "(RPKM)",
            cuttag = TRUE
        ),
        H4 = list(
            distUnit = 1,
            include = c("quantile"),
            log10 = FALSE,
            label = "H4 Cut&Tag",
            unit = "# Reads",
            enrichmentLabel = "H4 (log10 RPKM)",
            trackHeaderLabel = "H4 Cut&Tag",
            trackScoreLabel = "RPKM",
            class = "chromatin",
            valueLim = c(0,50),
            deltaLim = c(-3, 3),
            corrLim  = c(0, 2),
            corrUnit = "(RPKM)",
            cuttag = TRUE
        ),
        H4ac = list(
            distUnit = 1,
            include = c("quantile"),
            log10 = FALSE,
            label = "H4ac Cut&Tag",
            unit = "# Reads",
            enrichmentLabel = "H4ac (log10 RPKM)",
            trackHeaderLabel = "H4ac Cut&Tag",
            trackScoreLabel = "RPKM",
            class = "chromatin",
            valueLim = c(0,50),
            deltaLim = c(-3, 3),
            corrLim  = c(0, 2),
            corrUnit = "(RPKM)",
            cuttag = TRUE
        ),
        H3K27me3 = list(
            distUnit = 1,
            include = c("quantile"),
            log10 = FALSE,
            label = "H3K27me3 Cut&Tag",
            unit = "# Reads",
            enrichmentLabel = "H3K27me3 (log10 RPKM)",
            trackHeaderLabel = "H3K27me3 Cut&Tag",
            trackScoreLabel = "RPKM",
            class = "chromatin",
            valueLim = c(0,50),
            deltaLim = c(-3, 3),
            corrLim  = c(0, 2),
            corrUnit = "(RPKM)",
            cuttag = TRUE
        )
    )
)
scoreTypes$genome$gc_z <- scoreTypes$genome$gc
