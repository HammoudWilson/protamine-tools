
# set any constants/environment variables needed by the script
LOG10_M=100 # multiplier when rounding log10 RPKM values

# declare the input file(s) used to create the output BED file
# this BED file was created by nascent/bin and carries 1kb bins with Pro-seq RPKM values for the unstranded composite genome
#  chrom   start   end   cpm==rpkm
INPUT_FILE=${TASK_DIR}/atac_060225_v6.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# keep all bins as is with no merging, round log10 RPKM values
# this provides the most accurate local transcription level for a small query region at the expense of a large file
zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | # remove the header
awk '$1 ~ /-'${PRIMARY_GENOME}'/' | # filter to the primary genome and correct chrom names
sed 's/-'${PRIMARY_GENOME}'//' |
awk -v M=${LOG10_M} 'BEGIN{
    OFS = "\t";
    print "chrom", "start0", "end1", "name", "Log10_RPKM";
} {
    if($4 < 0.001) {
        log10_rpkm = -3;
    } else {
        log10_rpkm = int(log($4)/log(10) * M + 0.5) / M;
    }
    print $1, $2, $3, ".", log10_rpkm;
}'
