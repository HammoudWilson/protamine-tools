
# shared script for parsing active nascent transcription spans at different RPKM thresholds
# MIN_RPKM and TXED_FILTER are set in the calling script
LOG10_M=100 # multiplier when rounding log10 RPKM values

# declare the input file(s) used to create the output BED file
# this BED file was created by nascent/bin and carries 1kb bins with Pro-seq RPKM values for the unstranded composite genome
#  chrom   start   end   cpm==rpkm
INPUT_FILE=${TASK_DIR}/${DATA_NAME}.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# filter to actively transcribed bins, round log10 RPKM values, and merge to output BED file
# RPKM averaging is imperfect since adjacent txn units with different txn levels will be averaged together in one span
zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | # remove the header
awk '$1 ~ /-'${PRIMARY_GENOME}'/' | # filter to the primary genome and correct chrom names
sed 's/-'${PRIMARY_GENOME}'//' |
awk -v MIN_RPKM=${MIN_RPKM} 'BEGIN{
    OFS = "\t";
} {
    txed = $4 >= MIN_RPKM ? 1 : 0;
    print $1, $2, $3, txed, $4;
}' |
bedtools groupby -g 1,4, -c 2,3,5 -o min,max,mean |
awk -v M=${LOG10_M} -v TXED_FILTER=${TXED_FILTER} 'BEGIN{
    OFS = "\t";
    print "chrom", "start0", "end1", "name", "Log10_RPKM";
} $2 == TXED_FILTER {
    if($5 < 0.001) {
        log10_rpkm = -3;
    } else {
        log10_rpkm = int(log($5)/log(10) * M + 0.5) / M;
    }
    print $1, $3, $4, ".", log10_rpkm;
}'
