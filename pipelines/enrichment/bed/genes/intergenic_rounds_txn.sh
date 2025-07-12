# set any constants/environment variables needed by the script
LOG10_M=100

# declare the input file(s) used to create the output BED file
INPUT_FILE_1=${GENOMES_DIR}/${PRIMARY_GENOME}/annotations/mm39.gencode.M36.genes.bed.gz
INPUT_FILE_2=${TASK_DIR}/atac_060225_v6.nascent_transcriptome_unstranded.bed.gz

# name and excerpt the input file(s) to the job log
# be sure to include the redirection to STDERR (1>&2) on all log file lines!
echo "input file 1: ${INPUT_FILE_1}" 1>&2
zcat ${INPUT_FILE_1} | head 1>&2
echo 1>&2
echo "input file 2: ${INPUT_FILE_2}" 1>&2
zcat ${INPUT_FILE_2} | head 1>&2

# do the work
zcat ${INPUT_FILE_1} |
bedtools complement -i - -g ${PRIMARY_GENOME_FASTA}.fai | 
bedtools intersect -loj -a - -b <(
    zcat ${INPUT_FILE_2} |
    awk '$1 != "chrom"' | 
    awk '$1 ~ /mm39/' | 
    sed 's/-mm39//'
) |
awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, ($7 == "." ? 0 : $7) }' |
bedtools groupby -g 1,2,3 -c 4 -o mean | 
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
