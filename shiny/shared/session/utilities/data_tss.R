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

# add series data to a plot based on plot type
paTSS_add_series <- function(inserts, yOffset, ylim_series, seriesRange, config){
    switch(
        config$Plot_Inserts_As,
        endpoint_counts = {
            segments(
                x0 = inserts$starts$x - 1/4,
                x1 = inserts$starts$x - 1/4,
                y0 = yOffset,
                y1 = yOffset + pmin(inserts$starts$y, ylim_series[2]) / seriesRange,
                col = CONSTANTS$plotlyColors$blue,
                lwd = 1.5
            )
            segments(
                x0 = inserts$ends$x + 1/4,
                x1 = inserts$ends$x + 1/4,
                y0 = yOffset,
                y1 = yOffset + pmin(inserts$ends$y, ylim_series[2]) / seriesRange,
                col = CONSTANTS$plotlyColors$red,
                lwd = 1.5
            )
        },
        v_plot = {
            abline(h = yOffset + 146 / seriesRange, col = "grey")
            points(
                x = inserts$x, 
                y = inserts$y / seriesRange + yOffset, 
                col = CONSTANTS$plotlyColors$blue %>% addAlphaToColor(0.5), 
                pch = 16,
                cex = 0.75
            )
        },
        insert_spans = {
            abline(h = yOffset + 146 / seriesRange, col = "grey")
            segments(
                x0 = inserts$x1, 
                y0 = inserts$y / seriesRange + yOffset, 
                x1 = inserts$x2, 
                y1 = inserts$y / seriesRange + yOffset, 
                col = CONSTANTS$plotlyColors$blue %>% addAlphaToColor(0.5), 
                lwd = 1
            )
        }
    )
}

# convert ATAC inserts into XY values based on plot type
paTss_parse_inserts <- function(inserts, metadata, config, nInserts, wInserts){
    switch(
        config$Plot_Inserts_As,
        endpoint_counts = {
            bases <- data.table(x = config$coordStart1:config$coordEnd1)
            x <- list(
                starts = merge(
                    bases,
                    # inserts[, .(y = .N / nInserts), keyby = .(x = start0 + 1)],
                    inserts[, .(y = sum(startWeight) / wInserts), keyby = .(x = start0 + 1)],
                    all.x = TRUE,
                    all.y = FALSE,
                    sort = FALSE,
                    by = "x"
                ),
                ends = merge(
                    bases,
                    # inserts[, .(y = .N / nInserts), keyby = .(x = end1)],
                    inserts[, .(y = sum(endWeight) / wInserts), keyby = .(x = end1)],
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
        v_plot = {
            inserts[, .(
                x = start0 + 1 + (end1 - start0) / 2, # midpoint
                y = end1 - start0 - metadata$env$MIN_INSERT_SIZE + 1 # length of insert
            )]
        },
        insert_spans = {
            inserts[, .(
                x1 = start0 + 1,
                x2 = end1,
                y = end1 - start0 - metadata$env$MIN_INSERT_SIZE + 1 # length of insert
            )]
        }
    )
}

# load ATAC inserts in a browser window
getAtacInserts <- function(metadata, config, coord, bgzFileName){ # returns a list of sample-level score objects based on GC normalization
    insertsDir <- trimws(config$Inserts_Bgz_Dir)
    if(insertsDir == "") insertsDir <- metadata$env$INSERTS_BGZ_DIR
    bgzFile <- file.path(insertsDir, bgzFileName)
    req(file.exists(bgzFile))
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
    )
}
paTss_get_inserts <- function(metadata, coord, config){
    switch(
        config$Aggregate_By,
        sample = {
            sapply(metadata$samples$sample_name, function(sample_name) {
                footprint <- metadata$footprint$sample[[sample_name]]
                getAtacInserts(metadata, config, coord, footprint$bgzFileName) %>%
                paTss_parse_inserts(metadata, config, footprint$nInserts, footprint$wInserts)
            }, simplify = FALSE, USE.NAMES = TRUE)
        },
        stage = {
            sapply(metadata$stages, function(stage) {
                footprint <- metadata$footprint$stage[[stage]]
                getAtacInserts(metadata, config, coord, footprint$bgzFileName) %>%
                paTss_parse_inserts(metadata, config, footprint$nInserts, footprint$wInserts)
            }, simplify = FALSE, USE.NAMES = TRUE)
        },
        stage_type = {
            sapply(names(metadata$stageTypes), function(stageType) {
                footprint <- metadata$footprint$stageType[[stageType]]
                getAtacInserts(metadata, config, coord, footprint$bgzFileName) %>%
                paTss_parse_inserts(metadata, config, footprint$nInserts, footprint$wInserts)
            }, simplify = FALSE, USE.NAMES = TRUE)
        }
    )
}
