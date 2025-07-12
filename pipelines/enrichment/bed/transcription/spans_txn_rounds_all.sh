
# set any constants/environment variables needed by the script
LOG10_M=10 # multiplier when rounding log10 RPKM values; lower number here tries to merge more bins into fewer spans

# declare the input file(s) used to create the output BED file
# this BED file was created by nascent/bin and carries 1kb bins with Pro-seq RPKM values for the unstranded composite genome
#  chrom   start   end   cpm==rpkm
INPUT_FILE=${TASK_DIR}/atac_060225_v6.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# keep all bins, round log10 RPKM values, and merge to output BED file
# output spans do not correspond to any specific transcription unit, just spans of similar transcription levels
zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | # remove the header
awk '$1 ~ /-'${PRIMARY_GENOME}'/' | # filter to the primary genome and correct chrom names
sed 's/-'${PRIMARY_GENOME}'//' |
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
