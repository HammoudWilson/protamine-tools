# action:
#     import two ChIP bigwig files, one for each genome strand
#     as needed, liftover the input data to GENOME
#     reformat the data to create a genome-wide map that is:
#         _unstranded_, i.e., values reflect signal on both strands
#         binned to BIN_SIZE bp bins, for direct comparison to ATAC-seq data
# input:
#     two bigwig files
#     genome bins BED files
#     optional liftover chain file
# where, we assumme ChIP data are:
#     ???
# outputs:
#     ChIP signal map as a genome-wide bin BED file where score = cpm over both strands
#     where cpm is Counts Per Million reads

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(rtracklayer)
    library(data.table)
    library(parallel)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOME_BINS_BED',
        'DATA_FILE_PREFIX',
        'BIGWIG_FILE_FORWARD',
        'BIGWIG_FILE_REVERSE',
        'LIFTOVER_CHAIN',
        'CHIP_NAME'
    ),
    integer = c(
        'BIN_SIZE',
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2)
#=====================================================================================

#=====================================================================================
# do the work
#-------------------------------------------------------------------------------------
message("importing (and lifting over) bigwig files")
chip <- do.call(rbind, mclapply(c(env$BIGWIG_FILE_FORWARD, env$BIGWIG_FILE_REVERSE), function(bw_file) {
# nascent <- do.call(rbind, lapply(c(env$BIGWIG_FILE_FORWARD, env$BIGWIG_FILE_REVERSE), function(bw_file) {

    message(paste(" ", "importing", bw_file))
    bw_data <- import(bw_file, format = "bigWig") # a GRanges object
    if(env$LIFTOVER_CHAIN != "NA") {
        message(paste(" ", "lifting with ", env$LIFTOVER_CHAIN))
        chain <- import.chain(env$LIFTOVER_CHAIN)
        bw_data <- unlist(liftOver(bw_data, chain))
    }

    # convert to data.table
    # width is always 1, strand is always "*", since F and R files are provided
    dt <- as.data.table(bw_data) 
    rm(bw_data)
    dt[, .(chrom = as.character(seqnames), start0 = start - 1L, count = score)]
}, mc.cores = env$N_CPU))
# }))

message("creating unstranded ChIP signal map by counting ????")
chip[, start0 := floor(start0 / env$BIN_SIZE) * env$BIN_SIZE] # the start0 of the bin containing the read 3' end
chip <- chip[, .(count = sum(count, na.rm = TRUE)), keyby = .(chrom, start0)] # sum the counts per bin

message("merging into genome bins, including empty bins")
bins <- fread(env$GENOME_BINS_BED)
chip <- merge(
    bins[, .(chrom, start0, end1)], 
    chip, 
    by = c('chrom', 'start0'), 
    all.x = TRUE,
    sort = FALSE
)
chip[is.na(count), count := 0L]

message("converting counts to CPM")
chip[, cpm := count / sum(count, na.rm = TRUE) * 1e6]

message("writing chip output file")
chip <- chip[, .(chrom, start0, end1, cpm)]
fwrite(
    chip, 
    file = paste(env$DATA_FILE_PREFIX, env$CHIP_NAME, "chip_unstranded.bed.gz", sep = '.'),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = '\t'
)
message()
str(chip)
