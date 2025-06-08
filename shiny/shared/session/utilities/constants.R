# limit values
gcLimits <- c(0.25, 0.6) # optimal for mm39 genome; for normalization plots only, all other inherit from pipeline env

# empirically determine insert size values
meanNucleosomeFragSize   <- 205
meanDinucleosomeFragSize <- 375

# data types
refTypes <- c('genome', 'spike_in')
insertTypes <- c('all_inserts', 'intermediate')

# score types metadata
scoreTypes <- list(
    genome = list(
        # see definition of gc_z below, copied from gc
        gc = list(
            label = "GC",
            unit = "Percent",
            trackHeaderLabel = "GC Percent",
            trackSummaryLabel = "Z Score",
            class = "baseComposition",
            valueLim = gcLimits
        ),
        txn = list(
            label = "Round Spermatid PRO-seq",
            unit = "log10 CPM",
            trackHeaderLabel = "Round Spermatid PRO-seq",
            trackSummaryLabel = "log10 CPM",
            class = "transcription",
            valueLim = log10(c(1e-3, 1e3))
        )
    ),
    sample = list(
        gcrz_obs = list(
            distUnit = 0.1,
            include = c("quantile"),
            log10 = FALSE,
            label = "GC Residual",
            unit = "Z Score",
            enrichmentLabel = "GC Regression Observed (Z)",
            trackHeaderLabel = "GC Regression Observed",
            trackScoreLabel = "Z Score",
            class = "coverage",
            valueLim = c(-3, 3),
            deltaLim = c(-3, 3)
        ),
        gcrz_wgt = list(
            distUnit = 0.1,
            include = c("quantile"),
            log10 = FALSE,
            label = "GC Residual",
            unit = "Z Score",
            enrichmentLabel = "GC Regression Tn5 Weighted (Z)",
            trackHeaderLabel = "GC Regression Tn5 Weighted",
            trackScoreLabel = "Z Score",
            class = "coverage",
            valueLim = c(-3, 3),
            deltaLim = c(-3, 3)
        ),
        nrll = list(
            distUnit = 0.1,
            label = "Protamine Transition Enrichment",
            unit = "NRLL",
            enrichmentLabel = "Protamine Transition (NRLL)",
            trackHeaderLabel = "Protamine Transition Enrichment",
            trackScoreLabel = "NRLL",
            class = "insertSize",
            valueLim = c(-1.5, 1.5),
            deltaLim = c(-1.5, 0.5)
        )
    )
)
scoreTypes$genome$gc_z <- scoreTypes$genome$gc
