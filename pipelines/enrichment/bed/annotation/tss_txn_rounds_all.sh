
# set any constants/environment variables needed by the script
TSS_PADDING=100
LOG10_M=100 # multiplier when rounding log10 RPKM values

# declare the input file(s) used to create the output BED file
# this tss bed was created by nascent/tss and carries 1kb-padded gene TSS with Pro-seq RPKM values for the composite genome
#   chrom   start   end     name    gene_size       strand  upstream_rpkm   gene_start_rpkm tss_fold_increase
#   chr1-mm39       3142476 3144476 ENSMUSG00000102693.2/4933401J01Rik      1070    +       0.0730458263765274      0.106717927066958       1.46097227399225 
INPUT_FILE=${TASK_DIR}/atac_060225_v6.mm39_dm6.gencode.M36.tss.bed.gz

# name and excerpt the input file(s) to the job log
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# adjust to 100bp-padded TSS regions with log10 RPKM values, all TSS
# minimum log10 RPKM is -3.0 == log10(0.001)
zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | # remove the header
awk '$1 ~ /-'${PRIMARY_GENOME}'/' | # filter to the primary genome and correct chrom names
sed 's/-'${PRIMARY_GENOME}'//' |
awk -v PADDING=${TSS_PADDING} -v M=${LOG10_M} 'BEGIN{
    OFS = "\t";
    print "chrom", "start0", "end1", "name", "Log10_TSS_RPKM";
} {
    if($8 < 0.001) {
        log10_rpkm = -3;
    } else {
        log10_rpkm = int(log($8)/log(10) * M + 0.5) / M;
    }
    print $1, $2 + 1000 - PADDING, $3 - 1000 + PADDING, $4, log10_rpkm;
}'
