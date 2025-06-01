#----------------------------------------------------------------------
# server components for the modelParameters appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
modelParametersServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'modelParameters'
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
# data sources
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "dataSourceTable", 
    selection = "single"
) 
spermatidSamplesTable <- spermatidSamplesTableServer(
    "spermatidSamplesTable", 
    sourceId,
    selection = "single"
)
selectedSample <- reactive({
    spermatidSamplesTable$selectedSamples()
})
selectedSampleName <- reactive({
    selectedSample()$sample_name
})

#----------------------------------------------------------------------
# mappability plots
#----------------------------------------------------------------------
logit <- function(x, min = 0.0001, max = 0.9999) {
    x <- pmax(pmin(x, max), min) # clip to range
    log(x / (1 - x))
}
mappabilityPlot <- staticPlotBoxServer(
    "mappabilityPlot",
    maxHeight = "400px",
    # lines   = TRUE,
    # legend  = TRUE,
    # margins = TRUE,
    # title   = TRUE,
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        req(sourceId, selectedSampleName)
        startSpinner(session, message = "rendering mappability")
        I <- paCollateData_v5(sourceId)$bins[, sample(which(primary & included), 2000)]
        bin_n_raw <- paCollateData_v5(sourceId)$bin_n_raw[I, selectedSampleName]
        bin_mpp   <- paCollateData_v5(sourceId)$bin_mpp[  I, selectedSampleName]
        bin_mpp <- pmax(pmin(bin_mpp, 0.9999), 0.0001)
        # aggregate <- input$aggregateByStage
        mappabilityPlot$initializeFrame(
            xlim = logit(c(0.0001, 0.9999)),
            ylim = c(0, quantile(bin_n_raw, 0.975, na.rm = TRUE)),
            xlab = "Logit Bin Mappability",
            ylab = "Bin Read Count"
        )
        # abline(v = seq(100, 600, 100), col = "grey80") # nucleosome size boundaries
        # abline(v = MAPPABILITY_KMER_LENGTHS - 0.5, col = "grey80")
        mappabilityPlot$addPoints(
            x = logit(bin_mpp),
            y = bin_n_raw,
            cex = 0.25,
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
    }
)
mappabilityDistributionPlot <- staticPlotBoxServer(
    "mappabilityDistributionPlot",
    maxHeight = "400px",
    # lines   = TRUE,
    # legend  = TRUE,
    # margins = TRUE,
    # title   = TRUE,
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        req(sourceId, selectedSampleName)
        startSpinner(session, message = "rendering mappability distribution")
        I <- paCollateData_v5(sourceId)$bins[, which(primary & included)]
        bin_mpp <- round(logit(paCollateData_v5(sourceId)$bin_mpp[I, selectedSampleName]), 1)
        d <- data.table(bin_mpp = bin_mpp)[, .N, by = .(bin_mpp)][, f := N / sum(N, na.rm = TRUE)]
        mappabilityDistributionPlot$initializeFrame(
            xlim = logit(c(1/3, 1)),
            ylim = c(0, max(d$f) * 1.05),
            xlab = "Bin Mappability",
            ylab = "Frequency"
        )
        # abline(v = seq(100, 600, 100), col = "grey80") # nucleosome size boundaries
        # abline(v = MAPPABILITY_KMER_LENGTHS - 0.5, col = "grey80")
        mappabilityDistributionPlot$addPoints(
            x = d$bin_mpp,
            y = d$f,
            cex = 0.25,
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# GC plots
#----------------------------------------------------------------------
randomBinI <- reactive({
    sourceId <- sourceId()
    selectedSampleName <- selectedSampleName()
    req(sourceId, selectedSampleName)
    startSpinner(session, message = "setting random bins")
    bin_mpp <- paCollateData_v5(sourceId)$bin_mpp[, selectedSampleName]
    paCollateData_v5(sourceId)$bins[, sample(which(primary & included & bin_mpp >= 1/3), 2000)]
})
gcPlot <- staticPlotBoxServer(
    "gcPlot",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        req(sourceId, selectedSampleName)
        I <- randomBinI()
        startSpinner(session, message = "rendering GC plot")
        bin_n_raw <- paCollateData_v5(sourceId)$bin_n_raw[I, selectedSampleName]
        bin_gc    <- paCollateData_v5(sourceId)$bin_gc[   I, selectedSampleName]
        gcPlot$initializeFrame(
            xlim = gcLimits,
            ylim = c(0, quantile(bin_n_raw, 0.975, na.rm = TRUE)),
            xlab = "Bin Fraction GC",
            ylab = "Bin Read Count"
        )
        gcPlot$addPoints(
            x = bin_gc,
            y = bin_n_raw,
            cex = 0.25,
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
    }
)
gcPlot2 <- staticPlotBoxServer(
    "gcPlot2",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        req(sourceId, selectedSampleName)
        I <- randomBinI()
        bin_n_raw <- paCollateData_v5(sourceId)$bin_n_raw[I, selectedSampleName]
        bin_gc    <- paCollateData_v5(sourceId)$bin_gc[   I, selectedSampleName]
        bin_mpp   <- paCollateData_v5(sourceId)$bin_mpp[  I, selectedSampleName]
        gcPlot2$initializeFrame(
            xlim = gcLimits,
            ylim = c(0, quantile(bin_n_raw, 0.975, na.rm = TRUE)),
            xlab = "Bin Fraction GC",
            ylab = "Corr. Bin Read Count"
        )
        gcPlot2$addPoints(
            x = bin_gc,
            y = bin_n_raw / bin_mpp,
            cex = 0.25,
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# Tn5 plots
#----------------------------------------------------------------------
tn5Dist_all <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "parsing Tn5 sites")
    n_tn5_all <- rowSums(paCollateData_v5(sourceId)$smp_n_tn5)
    f_tn5_all <- n_tn5_all / sum(n_tn5_all, na.rm = TRUE)
    I <- order(f_tn5_all, decreasing = TRUE)
    cdf <- cumsum(f_tn5_all[I])
    I <- I[cdf <= 0.9999] # keep only the top 95% of Tn5 sites
    rank <- 1:length(I)
    data.table(
        rank = rank,
        kmer = paCollateData_v5(sourceId)$tn_5_kmers[I],
        f    = f_tn5_all[I],
        cdf  = cdf[rank],
        I    = I
    )
})
tn5Plot <- staticPlotBoxServer(
    "tn5Plot",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        tn5Dist_all <- tn5Dist_all()
        req(selectedSampleName, tn5Dist_all)
        startSpinner(session, message = "rendering Tn5 freq plot")
        smp_n_tn5 <- paCollateData_v5(sourceId)$smp_n_tn5[tn5Dist_all$I, selectedSampleName]
        smp_f_tn5 <- smp_n_tn5 / sum(smp_n_tn5, na.rm = TRUE)
        x_max <- nrow(tn5Dist_all)
        tn5Plot$initializeFrame(
            xlim = log10(c(1, x_max)),
            ylim = log10(c(min(tn5Dist_all$f[tn5Dist_all$f > 0]), max(smp_f_tn5, tn5Dist_all$f, na.rm = TRUE) * 1.05)),
            xlab = "Tn5 Kmer Rank",
            ylab = "Fraction of Tn5 Sites"
        )
        abline(h = 0, col = "grey80") # vertical line at 0
        abline(v = 1, col = "grey80") # vertical line at 0
        tn5Plot$addPoints(
            x = log10(1:x_max),
            y = log10(smp_f_tn5[1:x_max]),
            cex = 0.25,
            col = CONSTANTS$plotlyColors$blue %>% addAlphaToColor(0.1)
        )
        tn5Plot$addLines(
            x = log10(1:x_max),
            y = log10(tn5Dist_all$f[1:x_max]),
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
    }
)
tn5Diff_RS_ES <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "parsing RS vs. ES")
    n_tn5_all <- paCollateData_v5(sourceId)$smp_n_tn5
    n_tn5_RS <- rowSums(n_tn5_all[, 3:11])
    n_tn5_ES <- rowSums(n_tn5_all[, 15:24])
    f_tn5_RS <- n_tn5_RS / sum(n_tn5_RS, na.rm = TRUE)
    f_tn5_ES <- n_tn5_ES / sum(n_tn5_ES, na.rm = TRUE)
    I <- order(f_tn5_RS, decreasing = TRUE)
    cdf <- cumsum(f_tn5_RS[I])
    I <- I[cdf <= 0.99]
    rank <- 1:length(I)
    data.table(
        rank_RS = rank,
        kmer = paCollateData_v5(sourceId)$tn_5_kmers[I],
        x    = log10(f_tn5_RS[I]),
        y    = log10(f_tn5_ES[I]),
        f_RS = f_tn5_RS[I],
        f_ES = f_tn5_ES[I],
        I    = I
    )
})
tn5DiffPlot <- interactiveScatterplotServer(
    "tn5DiffPlot", 
    tn5Diff_RS_ES,
    accelerate = TRUE,
    mode = "markers",
    color = NA,
    pointSize = 1,
    xtitle = "log10 Frac. Tn5 Sites (RS)",
    ytitle = "log10 Frac. Tn5 Sites (ES)",
    unityLine = TRUE,
    selectable = "lasso"
)
expandKmers_simple <- function(kmers) {
    unname(sapply(kmers, function(kmer) {
        kmer <- strsplit(kmer, "")[[1]]
        paste(c(kmer [1:3], "nn", kmer [4:6], "nn", kmer [7:9]), collapse = "")
    }))
}
observeEvent(tn5DiffPlot$selected(), {
    selected <- tn5DiffPlot$selected()
    tn5Diff <- tn5Diff_RS_ES()
    # 'data.frame':   3829 obs. of  5 variables:
    # $ curveNumber: int  0 0 0 0 0 0 0 0 0 0 ...
    # $ pointNumber: int  1260 1610 1713 1730 1911 2057 2313 2857 3077 3306 ...      
    # $ x          : num  -4.33 -4.37 -4.38 -4.38 -4.4 ...
    # $ y          : num  -4.8 -4.79 -5.05 -5.04 -5.05 ...
    # $ customdata : logi  NA NA NA NA NA NA ...
    i <- selected$pointNumber + 1 # ploty pointNumber is 0-based
    fwrite(
        tn5Diff[i, .(
            rank_RS = rank_RS,
            kmer = expandKmers_simple(kmer),
            f_RS = f_RS,
            f_ES = f_ES
        )], 
        file = "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/hammoud/protamine/tn5Diff.csv", 
        append = FALSE, 
        quote = FALSE,
        sep = ",", 
        row.names = FALSE, 
        col.names = TRUE,
        scipen = getOption('scipen', 0L)
    )
})
tn5CdfPlot <- staticPlotBoxServer(
    "tn5CdfPlot",
    maxHeight = "400px",
    create = function() {
        sourceId <- sourceId()
        selectedSampleName <- selectedSampleName()
        tn5Dist_all <- tn5Dist_all()
        req(selectedSampleName, tn5Dist_all)
        startSpinner(session, message = "rendering Tn5 cdf plot")
        smp_n_tn5 <- paCollateData_v5(sourceId)$smp_n_tn5[, selectedSampleName]
        smp_f_tn5 <- smp_n_tn5 / sum(smp_n_tn5, na.rm = TRUE)
        cdf <- cumsum(smp_f_tn5[tn5Dist_all$I])
        x_max <- nrow(tn5Dist_all)
        tn5CdfPlot$initializeFrame(
            xlim = c(1, x_max),
            ylim = c(0, 1),
            xlab = "Tn5 Kmer Rank",
            ylab = "Cumulative Fraction of Tn5 Sites"
        )
        abline(h = 0:1, col = "grey80") # horizontal line at 0
        tn5CdfPlot$addPoints(
            x = 1:x_max,
            y = cdf,
            cex = 0.5,
            col = CONSTANTS$plotlyColors$blue %>% addAlphaToColor(0.25)
        )
        tn5CdfPlot$addLines(
            x = 1:x_max,
            y = tn5Dist_all$cdf,
            col = CONSTANTS$plotlyColors$black
        )
        stopSpinner(session)
    }
)
expandKmers <- function(kmers) {
    unname(sapply(kmers, function(kmer) {
        kmer_ <- chartr("ACGT", "TGCA", kmer)
        kmer  <- strsplit(kmer,  "")[[1]]
        kmer_ <- strsplit(kmer_, "")[[1]]
        paste(
            paste(c(kmer [1:2], "|", kmer [3], "nn", kmer [4:6], "nn", kmer [7], "-", kmer [8:9]), collapse = ""),
            paste(c(kmer_[1:2], "-", kmer_[3], "nn", kmer_[4:6], "nn", kmer_[7], "|", kmer_[8:9]), collapse = ""),
            sep = "<br>"
        )
    }))
}
tn5TableData <- reactive({
    tn5Dist_all <- tn5Dist_all()
    req(tn5Dist_all)
    stopSpinner(session)
    tn5Dist_all[1:100, .(
        rank = rank,
        kmer = kmer,
        site = expandKmers(kmer),
        f    = round(f,   5),
        cdf  = round(cdf, 5)
    )]
})
bufferedTableServer(
    "tn5Table",
    id,
    input,  
    tn5TableData,
    selection = 'none'
)

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
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
