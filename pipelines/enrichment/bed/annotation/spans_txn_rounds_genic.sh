# set any constants/environment variables needed by the script
LOG10_M=100 # multiplier when rounding log10 RPKM values

# declare the input file(s) used to create the output BED file
# the first BED6 file is the gene-level annotation file for the mm39 primary genome created by genome/download
INPUT_FILE_1=${GENOMES_DIR}/${PRIMARY_GENOME}/annotations/${PRIMARY_GENOME}.gencode.M36.genes.bed.gz
# the second BED file was created by nascent/bin and carries 1kb bins with Pro-seq RPKM values for the unstranded composite genome
#  chrom   start   end   cpm==rpkm
INPUT_FILE_2=${TASK_DIR}/${DATA_NAME}.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
# be sure to include the redirection to STDERR (1>&2) on all log file lines!
echo "input file 1: ${INPUT_FILE_1}" 1>&2
zcat ${INPUT_FILE_1} | head 1>&2
echo 1>&2
echo "input file 2: ${INPUT_FILE_2}" 1>&2
zcat ${INPUT_FILE_2} | head 1>&2

# intersect genes with bins and average RPKM to get an approximate transcription level for each gene
zcat ${INPUT_FILE_1} |
bedtools intersect -loj -a - -b <(
    zcat ${INPUT_FILE_2} |
    awk '$1 != "chrom"' | 
    awk '$1 ~ /-'${PRIMARY_GENOME}'/' | 
    sed 's/-'${PRIMARY_GENOME}'//'
) |
awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4, ($10 == "." ? 0 : $10) }' | 
bedtools groupby -g 1,2,3,4 -c 5 -o mean | 
awk -v M=${LOG10_M} 'BEGIN{
    OFS = "\t";
    print "chrom", "start0", "end1", "name", "Log10_Gene_RPKM";
} {
    if($5 == "." || $5 < 0.001) {
        log10_rpkm = -3;
    } else {
        log10_rpkm = int(log($5)/log(10) * M + 0.5) / M;
    }
    print $1, $2, $3, $4, log10_rpkm;
}'
