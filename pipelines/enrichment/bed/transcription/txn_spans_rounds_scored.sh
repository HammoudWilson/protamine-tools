
# set any constants/environment variables needed by the script
export LOG10_M=10 # multiplier when rounding log10 RPKM values

# declare the input file(s) used to create the output BED file
export INPUT_FILE=${TASK_DIR}/atac_060225_v6.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# do whatever work is needed to create the output BED file as ${REGIONS_TYPE_BED} from the input file(s)
# R, bedtools, samtools, etc. are available for use if needed
# output file REGIONS_TYPE_BED must:
#     be written to STDOUT (don't save to file)
#     have a header line of column names
#     be proper BED format with at least 4 columns:
#         1. chromosome name,             column name 'chrom', e.g., 'chr3'
#         2. 0-referenced start position, column name 'start' or 'start0'
#         3. 1-referenced end position,   column name 'end'   or 'end1'
#         4. name of the region,          column name 'name', use . if name not applicable
#         optionally, BED files may have a 5th score column to use for score-type enrichment analysis
#             the name you give the to score column will be used as the axis label for enrichment plots
#             BED files without a score column can still be used for fraction-overlap-type enrichment analysis
#         you may include additional columns after the 5th score column, but they will be ignored during enrichment analysis

zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | 
awk '$1 ~ /mm39/' | 
sed 's/-mm39//' |
awk -v M=${LOG10_M} 'BEGIN{
    OFS = "\t";
} {
    if($4 < 0.001) {
        log10_rpkm = -3;
    } else {
        log10_rpkm = int(log($4)/log(10) * M + 0.5) / M;
    }
    print $1, $2, $3, log10_rpkm;
}' |
bedtools groupby -g 1,4, -c 2,3 -o min,max |
awk 'BEGIN{
    OFS = "\t";
    print "chrom", "start0", "end1", "name", "Log10_RPKM";
} {
    print $1, $3, $4, ".", $2;
}'
