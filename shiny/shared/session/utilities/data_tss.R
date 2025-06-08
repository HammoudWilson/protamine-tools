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

# add one group's data to a footprint track plot based on plot type
paTSS_add_series <- function(inserts, yOffset, ylim_series, seriesRange, seriesName, metadata, config){
    sizeH <- c(metadata$env$MIN_INSERT_SIZE, meanNucleosomeFragSize, meanDinucleosomeFragSize) # horizontal lines for size-based traces
    sizeH_ <- yOffset + sizeH / seriesRange # as above now in virtual axis units
    sizedColor <- function() {
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
        }
    )
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
