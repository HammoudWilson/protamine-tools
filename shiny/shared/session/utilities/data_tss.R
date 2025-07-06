# load and format TSS data
paTss_frags <- function(sourceId){
    startSpinner(session, message = "loading TSS inserts")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "tssFrags", 
        ttl = CONSTANTS$ttl$month
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# load footprint metadata
paTss_footprint <- function(sourceId){
    startSpinner(session, message = "loading TSS footprint")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "footprint", 
        ttl = CONSTANTS$ttl$month
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}

# load ab_initio data
paTss_ab_initio <- function(sourceId){
    startSpinner(session, message = "loading nuc chains")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "ab_initio", 
        ttl = CONSTANTS$ttl$month
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}
paTss_coordCols <- c("chrom", "start0", "end1")
paTss_clusterCols <- c(paTss_coordCols, "k_scaled", "k_unscaled")
paTss_appRPKMCols <- c("stage_mean", "min_RPKM", "max_RPKM", "delta_RPKM")
paTss_profileCols <- c("mean_stage", "max_stage", "quantile", paTss_appRPKMCols)
paTss_profile_ab_initio <- function(regions, stages){
    nStages <- length(stages)
    regions$min_RPKM   <- apply(regions[, ..stages], 1, min)
    regions$max_RPKM   <- apply(regions[, ..stages], 1, max)
    regions$delta_RPKM <- (regions$max_RPKM - regions$min_RPKM) %>% round(2)
    regions$stage_mean <- apply(regions[, ..stages], 1, function(rpkm) weighted.mean(1:nStages, rpkm))
    regions$max_stage  <- stages[apply(regions[, ..stages], 1, which.max)]
    regions$mean_stage <- stages[round(regions$stage_mean)]
    regions$quantile <- as.integer(cut(
        regions$stage_mean, 
        breaks = quantile(regions$stage_mean, probs = seq(0, 1, 0.05)), 
        include.lowest = TRUE
    ))
    regions$stage_mean <- regions$stage_mean %>% round(2)
    regions
}
paTss_dinuc_regions <- function(sourceId){
    startSpinner(session, message = "loading dinuc regions")
    x <- protaminerCache$get(
        "paTss_dinuc_regions",
        key = sourceId,
        create = "asNeeded",
        createFn = function(...) {
            startSpinner(session, message = "processing dinuc regions")
            fp <- paTss_footprint(sourceId)
            ai <- paTss_ab_initio(sourceId)
            rpkmCols <- paste(ai$stages, "rpkm", sep = '_')
            coordCols <- c("chrom", "start0", "end1")
            regions <- ai$regions[index_stage == "overlap_group", .SD, .SDcols = c(paTss_clusterCols, rpkmCols)]
            setnames(regions, c(paTss_clusterCols, ai$stages))
            regions <- paTss_profile_ab_initio(regions, ai$stages)
            regions[, .SD, .SDcols = c(paTss_clusterCols, paTss_profileCols, ai$stages)]
        }
    )
    stopSpinner(session)
    x$value
}

# get collated lists of CTCF motifs for the reference genome
ctcf_motif_sources <- list(
    # AH104745 | mm39.JASPAR2022_CORE_vertebrates_non_redundant_v2.RData == CTCF sites detected using the all three CTCF PWMs
    # AH104747 | mm39.MA0139.1.RData  == CTCF sites detected using the most popular MA0139.1 CTCF PWM
    mm39 = list(
        dataprovider = "JASPAR 2022",
        database     = "AH104747" 
    )
)
paTss_ctcf_motifs <- function(reference){
    ctcf_motif_source <- ctcf_motif_sources[[reference$genome$genome]]
    req(ctcf_motif_source)
    startSpinner(session, message = "loading CTCF motifs")
    x <- protaminerCache$get(
        "paTss_ctcf_motifs",
        keyObject = ctcf_motif_source,
        create = "asNeeded",
        createFn = function(...) {
            # https://dozmorovlab.github.io/CTCF/
            ctcf <- AnnotationHub::subset(
                AnnotationHub::AnnotationHub(), 
                preparerclass == "CTCF" & 
                species == reference$genome$scientificName & 
                genome  == reference$genome$genome & 
                dataprovider == ctcf_motif_source$dataprovider
            )
            ctcf <- ctcf[[ctcf_motif_source$database]]
            ctcf <- data.table(
                chrom = as.character(seqnames(ctcf)),
                start0 = start(ctcf) - 1,
                end1 = end(ctcf),
                strand = as.character(strand(ctcf))
            )[, .(
                chrom,
                pos1 = start0 + (end1 - start0) / 2 + 1, # motifs are variable lengths e.g., 34 35 19 34 35 35 19,
                strand # CTCF motifs are stranded
            )]
            setkey(ctcf, pos1)
            ctcf
        }
    )
    stopSpinner(session)
    x$value
}

# add one group's data to a footprint track plot based on plot type
paTSS_add_series <- function(inserts, yOffset, ylim_series, seriesRange, seriesName, metadata, config){
    sizeH <- c(metadata$env$MIN_INSERT_SIZE, meanNucleosomeFragSize, meanDinucleosomeFragSize) # horizontal lines for size-based traces
    sizeH_ <- yOffset + sizeH / seriesRange # as above now in virtual axis units
    sizedColor <- function(transparent = TRUE) {
        color <- if(config$Color_Inserts_By == "stage"){
            switch(
                config$Aggregate_By,
                sample     = getSampleColorsByStage(metadata$samples, metadata$samples[sample_name == seriesName]),
                stage      = getStageColor(metadata, seriesName),
                stage_type = getStageTypeColor(metadata, seriesName)
            )
        } else { 
            CONSTANTS$plotlyColors[[config$Color_Inserts_By]] 
        } 
        if(!transparent) return(color)
        color %>% addAlphaToColor(config$Insert_Color_Alpha)
    }
    switch(
        config$Plot_Inserts_As,
        
        # insert start counts/weights per base are plotted with a positive direction from zero in blue
        # insert end   counts/weights per base are plotted with a negative direction from zero in purple
        endpoint_counts = {
            zeroLine <- yOffset + 0.5 # create a zero axis in the middle of this groups plotting zone
            abline(h = zeroLine, col = "grey")
            maxY <- ylim_series[2] * 0.95 # make sure group have a bit of white space between them
            axis(2, at = zeroLine + c(0.25, -0.25), labels = c("starts", "ends"), tick = FALSE, las = 1)
            segments(
                x0 = inserts$starts$x,
                x1 = inserts$starts$x,
                y0 = zeroLine,
                y1 = zeroLine + (pmin(inserts$starts$y, maxY) / seriesRange) / 2,
                col = CONSTANTS$plotlyColors$blue,
                lwd = config$Endpoint_Line_Width
            )
            segments(
                x0 = inserts$ends$x,
                x1 = inserts$ends$x,
                y0 = zeroLine,
                y1 = zeroLine - (pmin(inserts$ends$y, maxY) / seriesRange) / 2,
                col = CONSTANTS$plotlyColors$purple,
                lwd = config$Endpoint_Line_Width
            )
        },

        # V plot and insert span plots show the insert size in a positive direction from zero
        # maintain axis to zero with a horizontal line at MIN_INSERT_SIZE; important to interpreting patterns
        # allow user to choose the color pattern, point/line size and opacity
        v_plot = {
            abline(h = sizeH_, col = "grey")
            axis(2, at = sizeH_, labels = sizeH, tick = FALSE, las = 1)
            points(
                x = inserts$x,
                y = inserts$y / seriesRange + yOffset,
                col = sizedColor(),
                pch = 16,
                cex = config$V_Plot_Point_Size
            )
        },
        insert_spans = {
            abline(h = sizeH_, col = "grey")
            axis(2, at = sizeH_, labels = sizeH, tick = FALSE, las = 1)
            segments(
                x0 = inserts$x1, 
                y0 = inserts$y / seriesRange + yOffset, 
                x1 = inserts$x2, 
                y1 = inserts$y / seriesRange + yOffset, 
                col = sizedColor(), 
                lwd = config$Insert_Line_Width
            )
        },

        # nucleosome bias is a simple XY plot
        nucleosome_bias = {
            abline(h = yOffset + 0.5, col = "grey")
            axis(2, at = (c(-0.9, 0, 0.9) + 1) / 2 + yOffset, labels = c("-1", "0", "1"), tick = FALSE, las = 1)
            # points(
            #     x = inserts$x,
            #     y = inserts$y / seriesRange + yOffset,
            #     col = sizedColor(FALSE),
            #     pch = 16,
            #     cex = 1
            # )
            lines(
                x = inserts$x,
                y = inserts$y / seriesRange + yOffset,
                col = sizedColor(FALSE),
                lty = 1,
                lwd = 1.5
            )
            # fd <- (inserts$y - 1) / 0.9
            # window_fn <- get(config$FFT_Window_Function, envir = loadNamespace("signal"))
            # fd <- fd * window_fn(length(fd)) # apply a window to the bias signal
            # fd <- abs(Re(fft(fd)))
            # # fd <- 1 / length(fd) * fd**2
            # fd <- fd / max(fd, na.rm = TRUE) # normalize the Fourier transform to [0, 1]
            # x <- config$coordStart1 + seq(0, length(fd) - 1) * config$Bias_Step_Size
            # segments(
            #     x0 = x, # * config$Bias_Step_Size,
            #     x1 = x,
            #     y0 = yOffset,
            #     y1 = fd + yOffset,
            #     col = sizedColor(FALSE),
            #     lwd = 1
            # )
        }
    )
}

# create an inserts plot, either in browser track or in a static plot
paTss_add_inserts_plot <- function(
    sourceId, coord, metadata, config, 
    nSeries, seriesNames, ylim_series, seriesRange, 
    inserts
){
    # overplot the called nucleosome chains as rectangles
    if(config$Show_Nucleosome_Chains && config$Aggregate_By == "stage"){
        regions <- paTss_ab_initio(sourceId)$regions[
            chrom  == coord$chromosome & 
            start0 <  coord$end & # wider than the plotted spans, includes the analysis flanks
            end1   >= coord$start
        ]
        if(nrow(regions) > 0){

            # overplot the called overlap groups across all stages
            regions_group <- regions[index_stage == "overlap_group"]
            nGroups <- nrow(regions_group)
            if(nGroups > 0) rect(
                regions_group$start0 + 1, 
                rep(0, nGroups),
                regions_group$end1, 
                rep(nSeries, nGroups), 
                border = CONSTANTS$plotlyColors$grey,
                lwd = config$Overlay_Line_Width
            )

            # overplot the called nucleosome chains by stage
            for(stageI in 1:nSeries){
                yOffset <- nSeries - stageI
                regions_stage <- regions[index_stage == seriesNames[stageI]]
                nChains <- nrow(regions_stage)
                if(nChains == 0) next
                rect(
                    regions_stage$start0 + 1, 
                    rep(yOffset + 0.025, nChains),
                    regions_stage$end1, 
                    rep(yOffset + 0.975, nChains), 
                    border = CONSTANTS$plotlyColors$red,
                    lwd = config$Overlay_Line_Width
                )
                for(chainI in 1:nChains){
                    nuc_starts0 <- as.integer(unlist(strsplit(regions_stage$nuc_starts0[chainI], ",")))
                    nNucs <- length(nuc_starts0)
                    rect(
                        nuc_starts0 + 1, 
                        rep(yOffset + 0.025, nNucs),
                        nuc_starts0 + 147, 
                        rep(yOffset + 0.975, nNucs),
                        border = NA, #CONSTANTS$plotlyColors$red,
                        col = CONSTANTS$plotlyColors$red %>% addAlphaToColor(0.1),
                        lwd = config$Overlay_Line_Width + 0.5
                    )
                }
            }
        }
    }

    # plot the inserts
    for(i in 1:nSeries){ 
        yOffset <- (nSeries - i)
        abline(h = yOffset, col = "black")
        paTSS_add_series(inserts[[i]], yOffset, ylim_series, seriesRange, seriesNames[i], metadata, config)
        text(coord$start, yOffset + 0.8, seriesNames[i], pos = 4, cex = 1.25)
    }
}

# convert ATAC inserts into XY values based on plot type for a specific sample/stage/stageType group 
paTss_parse_inserts <- function(inserts, metadata, config, footprint){
    switch(
        config$Plot_Inserts_As,

        # counts tracks show normalized weights per base for insert starts and ends independently
        # the user can choose to plot the counts as Tn5-preference-adjusted weights or as observed counts
        # normalization can be determined 
        #   per window for each group, so all groups will have some high values
        #   over the entire genome, so some groups could be high in a window while others are low
        # the y-axis is then scaled to the values based on quantile calculated over _all_ groups (see track ylim)
        endpoint_counts = {
            isWeights    <- config$Endpoint_Counts_As    == "weights"
            normToWindow <- config$Endpoint_Normalize_To == "window"
            bases <- data.table(x = config$coordStart1:config$coordEnd1) # ensure that all bases have at least zero values
            x <- list(
                starts = merge(
                    bases,
                    # footprint$xInserts was calculated over entire genome by pipeline
                    # sum(xWeight) or insert count calculated per window here (the only data we have access to)
                    if(isWeights){ 
                        wInserts <- if(normToWindow) inserts[, sum(startWeight)] else footprint$wInserts
                        inserts[, .(y = sum(startWeight) / wInserts), keyby = .(x = start0 + 1)]
                    } else {
                        nInserts <- if(normToWindow) inserts[, .N] else footprint$nInserts
                        inserts[, .(y = .N / nInserts), keyby = .(x = start0 + 1)]
                    },
                    all.x = TRUE,
                    all.y = FALSE,
                    sort = FALSE,
                    by = "x"
                ),
                ends = merge(
                    bases,
                    if(isWeights){
                        wInserts <- if(normToWindow) inserts[, sum(endWeight)] else footprint$wInserts
                        inserts[, .(y = sum(endWeight) / wInserts), keyby = .(x = end1)]
                    } else {
                        nInserts <- if(normToWindow) inserts[, .N] else footprint$nInserts
                        inserts[, .(y = .N / nInserts), keyby = .(x = end1)]
                    },
                    all.x = TRUE,
                    all.y = FALSE,
                    sort = FALSE,
                    by = "x"
                )
            )
            x$starts[, y := ifelse(is.na(y), 0, y)]
            x$ends[,   y := ifelse(is.na(y), 0, y)]
            x
        },

        # V plot and insert span tracks use simple conversion of inserts into plottable XY values per insert
        # while enforcing the maximum insert size expected by the plot to avoid group overruns
        # these plots show all inserts and have no easy way to normalize the values between groups
        v_plot = {
            inserts[end1 - start0 <= config$Max_Insert_Size, .(
                x = start0 + 1 + (end1 - start0) / 2, # midpoint
                y = end1 - start0 # length of insert
            )]
        },
        insert_spans = {
            inserts[end1 - start0 <= config$Max_Insert_Size, .(
                x1 = start0 + 1,
                x2 = end1,
                y = end1 - start0
            )]
        },

        # nucleosome bias plot the relative signal for mononucleosomes vs. intervening nucleosome-free regions
        nucleosome_bias = {
            dt <- inserts[, {
                size    <- end1 - start0
                start1  <- start0 + 1
                center1 <- start1 + size / 2
                is_mono_nuc <- size > 147 & size <= 320
                is_di_nuc   <- size > 320 & size <= 485
                is_sub_nuc  <- size < 147
                if(config$Bias_Method == "centers_only" && config$Centers_Only_Weighting == "Tn5_weight"){
                    weight <- pmax(startWeight, 5) + pmax(endWeight, 5)
                }
                m <- sapply(seq(config$coordStart1, config$coordEnd1 - config$Bias_Window_Size, by = config$Bias_Step_Size), function(windowStart1){
                    windowEnd1 <- windowStart1 + config$Bias_Window_Size - 1
                    center_is_in_window <- between(center1, windowStart1, windowEnd1)
                    switch(
                        config$Bias_Method,
                        centers_and_endpoints = {
                            start_is_in_window  <- between(start1,  windowStart1, windowEnd1)
                            end_is_in_window    <- between(end1,    windowStart1, windowEnd1)
                            nuc <- sum(is_mono_nuc & center_is_in_window) * 2
                            nfr <- sum(is_di_nuc   & center_is_in_window) * 2 + sum(start_is_in_window) + sum(end_is_in_window)
                        },
                        centers_only = {
                            switch(
                                config$Centers_Only_Weighting,
                                count = {
                                    nuc <- sum(is_mono_nuc & center_is_in_window)
                                    nfr <- sum((is_sub_nuc | is_di_nuc) & center_is_in_window)
                                },
                                Tn5_weight = {
                                    nuc <- ifelse( is_mono_nuc             & center_is_in_window, weight, 0) %>% sum()
                                    nfr <- ifelse((is_sub_nuc | is_di_nuc) & center_is_in_window, weight, 0) %>% sum()
                                }
                            )
                        }
                    )
                    denom <- nuc + nfr
                    bias <- if(denom == 0) NA else (nuc - nfr) / denom
                    c(
                        x = windowStart1 + config$Bias_Window_Size / 2, # midpoint of the window
                        y = bias
                    )
                })
                if(config$Center_Bias_Plots) m[2, ] <- m[2, ] - median(m[2, ], na.rm = TRUE) # center the bias around zero
                if(config$Scale_Bias_Plots)  m[2, ] <- m[2, ] / max(abs(m[2, ]), na.rm = TRUE) # scale the bias to [-1, 1]
                .(
                    x = m[1, ],
                    y = m[2, ] * 0.9 + 1
                )
            }]
        }
    )
}

# collect ATAC inserts in a browser window for a specific sample/stage/stageType group 
getAtacInserts <- function(metadata, config, coord, footprint){ # returns a list of sample-level score objects based on GC normalization
    
    # determine the bgz file from which to load inserts for this grouped footprint as created by the pipeline
    insertsDir <- trimws(config$Inserts_Bgz_Dir)
    if(insertsDir == "") insertsDir <- metadata$env$INSERTS_BGZ_DIR
    bgzFile <- file.path(insertsDir, footprint$bgzFileName)
    req(file.exists(bgzFile))

    # use tabix to load the inserts in the requested coordinate range from the bgz file
    tabix <- getCachedTabix(bgzFile)
    getTabixRangeData(
        tabix, 
        coord, 
        col.names = c(
            "chrom", 
            "start0", 
            "end1", 
            "startWeight",
            "endWeight"
        ), 
        colClasses = c(
            "character",
            "integer",
            "integer", 
            "numeric",
            "numeric"
        )
    ) %>% 

    # pass the call to parse the inserts by plot type
    paTss_parse_inserts(metadata, config, footprint)
}

# load ATAC inserts in a browser window, grouped by sample/stage/stageType as requested
paTss_get_inserts <- function(metadata, coord, config){
    switch(
        config$Aggregate_By,
        sample = {
            sapply(metadata$samples$sample_name, function(sample_name) {
                footprint <- metadata$footprint$sample[[sample_name]]
                getAtacInserts(metadata, config, coord, footprint)
            }, simplify = FALSE, USE.NAMES = TRUE)
        },
        stage = {
            sapply(metadata$stages, function(stage) {
                footprint <- metadata$footprint$stage[[stage]]
                getAtacInserts(metadata, config, coord, footprint)
            }, simplify = FALSE, USE.NAMES = TRUE)
        },
        stage_type = {
            sapply(names(metadata$stageTypes), function(stageType) {
                footprint <- metadata$footprint$stageType[[stageType]]
                getAtacInserts(metadata, config, coord, footprint)
            }, simplify = FALSE, USE.NAMES = TRUE)
        }
    )
}
