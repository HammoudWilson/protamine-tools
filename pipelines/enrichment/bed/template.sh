
# set any constants/environment variables needed by the script
ABC=123

# declare the input file(s) used to create the output BED file
# use comments here to describe the input file(s) source and content
INPUT_FILE=/path/to/input/file

# name and excerpt the input file(s) to the job log
# be sure to include the redirection to STDERR (1>&2) on all log file lines!
# choose `head` or `zcat` depending on whether the input file is compressed
echo "input file: ${INPUT_FILE}" 1>&2
head ${INPUT_FILE} 1>&2 
# zcat ${INPUT_FILE} | head 1>&2

# do whatever work is needed to create the output BED file as ${REGIONS_TYPE_BED} from the input file(s)
# R, bedtools, samtools, etc. are available for use if needed
# output file REGIONS_TYPE_BED must:
#     be written to STDOUT (don't save to file)
#     have a header line of column names
#     be proper BED format with at least 4 columns:
#         1. chromosome name, e.g., 'chr3', column name 'chrom'
#         2. 0-referenced start position,   column name 'start' or 'start0'
#         3. 1-referenced end position,     column name 'end'   or 'end1'
#         4. name of the region,            column name 'name'
#            use '.' in the name field if region names are not applicable
#     optionally, BED files may have a 5th score column to use for score-type enrichment analysis
#       the name you give the to score column will be used as the axis label for enrichment plots
#       BED files without a score column can still be used for fraction-overlapping enrichment analysis
#       if you include additional columns after the 5th score column, they will be dropped before writing to file
# note that overlapping spans in the BED stream will be merged into one span with the max score before writing to file
# a consequence of everything above is that BED spans/regions for ATAC-seq analysis are unstranded
