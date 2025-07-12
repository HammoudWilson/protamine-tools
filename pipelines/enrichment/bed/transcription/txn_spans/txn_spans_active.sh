
# shared script for parsing active nascent transcription spans at different RPKM thresholds
export LOG10_M=100 # multiplier when rounding log10 RPKM values

# declare the input file(s) used to create the output BED file
INPUT_FILE=${TASK_DIR}/atac_060225_v6.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# do whatever work is needed to create the output BED file as ${REGIONS_TYPE_BED} from the input file(s)
zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | 
awk '$1 ~ /mm39/' | 
sed 's/-mm39//' |
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
    log10_rpkm = int(log($5)/log(10) * M + 0.5) / M;
    print $1, $3, $4, ".", log10_rpkm;
}'
