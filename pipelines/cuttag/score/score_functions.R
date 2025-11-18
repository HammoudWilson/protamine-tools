# utilities for parsing bin scores and distributions
# overrwrites only the objects and functions that change relative to ATAC-seq scores

# score types metadata; minimal information only as required to calculate score distributions
scoreTypes <- list(
    sample = list(
        H2B = list(
            name = "H2B",
            distUnit = 1,
            include = c("quantile"), 
            log10 = FALSE,
            minValue = 1e-3
        ),
        H4 = list(
            name = "H4",
            distUnit = 1,
            include = c("quantile"), 
            log10 = FALSE,
            minValue = 1e-3
        ),
        H3K27me3 = list(
            name = "H3K27me3",
            distUnit = 1,
            include = c("quantile"), 
            log10 = FALSE,
            minValue = 1e-3
        ),
        H4ac = list(
            name = "H4ac",
            distUnit = 1,
            include = c("quantile"), 
            log10 = FALSE,
            minValue = 1e-3
        )
    )
)

# sample-level score type functions
# must return a vector the length of all composite bins, but scores only needed for primary genome bins
# or return NA if score type is not applicable to the sample
get_cuttag <- function(sample_name_){

# /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/Modules/biotoolbox/bin/bam2wig.pl 
# --out 26615R_histone/day56_ES_H2b.count.bw 
# --qual 10 --nosecondary --noduplicate --nosupplementary 
# --cpu 12 --bw --bwapp /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/Modules/biotoolbox/bin/wigToBigWig 
# --pe --mid --minsize 40 --maxsize 400 --chrskip 'chrM|chrUn|random' --rpm --scale 16.867858  
# --in 26615R_histone/day56_ES_H2b.dedup.bam  2>&1 > 26615R_histone/day56_ES_H2b.count.bam2wig.out.txt

    sample <- cuttagSamples[sample_name == sample_name_]
    m_reads <- sample$n_reads / 1e6
    bw_file <- file.path(env$INPUT_DIR, sample$sub_directory, "Count", sample$bigwig_file)
    cuttag <- {
        bw_data <- import(bw_file, format = "bigWig") # a GRanges object
        # convert to data.table
        # width is always 1, strand is always "*", since analysis is unstranded
        # chrom names in bigwig files expected to match PRIMARY_GENOME with chrom names like chr1, chrX
        dt <- as.data.table(bw_data) 
        rm(bw_data)
        # undo the count scaling in the bam2wig command to yield read count per called/non-zero base
        # as expected, most bases have a count of zero (not present in bigwig) or one read
        # depending on the antibody target and read depth, the majority of bases might have a count of zero
        dt[, .(chrom = as.character(seqnames), start0 = start - 1L, n_reads = as.integer(round(score / sample$scale_factor * m_reads)))]
    }
    cuttag[, start0 := floor(start0 / env$BIN_SIZE) * env$BIN_SIZE] # the start0 of the bin containing the read midpoint
    cuttag <- cuttag[, .(n_reads = sum(n_reads, na.rm = TRUE)), keyby = .(chrom, start0)] # sum the counts per bin across all bp
    cuttag[, chrom := paste0(chrom, "-", env$PRIMARY_GENOME)]
    cuttag <- merge(
        bins[, .(chrom, start0)], 
        cuttag, 
        by = c('chrom', 'start0'), 
        all.x = TRUE,
        sort = FALSE
    )
    cuttag[is.na(n_reads), n_reads := 0]

    # must return a vector the length of all composite bins, but scores only needed for primary genome bins
    # as with txn, cpm == rpkm for 1kb bins
    # as above, even at the bin level a majority of bins might have a count of zero
    cuttag$n_reads
}

# analyze and aggregate distributions of different bin scores
# all scores are expected to be on a comparable scale between samples
analyzeSampleScores <- function(scoreType, antibody_target_){
    sample_names <- cuttagSamples[antibody_target == antibody_target_, sample_name]
    x <- mclapply(sample_names, function(sample_name){
    # x <- lapply(sample_names, function(sample_name){
        message(paste("   ", "analyzeSampleScores", sample_name))
        analyzeScoreDist(get_cuttag(sample_name), scoreType, sample_name, return_scores = TRUE)
    }, mc.cores = env$N_CPU)
    # })
    names(x) <- sample_names
    x
}
aggregateSampleScores <- function(sampleScores, scoreType, antibody_target_){

    # aggregate scores by spermatid stage
    allStages <- cuttagSamples[antibody_target == antibody_target_, unique(stage)]
    by_stage <- mclapply(allStages, function(stage_){
    # by_stage <- lapply(allStages, function(stage_){
        sample_names <- cuttagSamples[antibody_target == antibody_target_ & stage == stage_, sample_name]
        if(length(sample_names) == 1){
            message(paste("   ", "single-sample stage", stage_))
            x <-sampleScores[[sample_names]]
            x$scores <- NULL
            x
        } else {
            message(paste("   ", "aggregateSampleScores by_stage", stage_))
            aggregateAndAnalyzeScores(sampleScores, sample_names, scoreType, stage_, rowAggFn = rowSums)
        }
    }, mc.cores = env$N_CPU)
    # })
    names(by_stage) <- allStages
    list(
        by_stage     = by_stage
    )
}
