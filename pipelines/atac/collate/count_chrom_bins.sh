# action:
#     count reads in genome bins for a single chrom in a single bam file
# input:
#     $BAM_FILE     as arg1
#     $CHROM        as arg2
#     $GENOME_FASTA_SHM, on /dev/shm
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $BIN_SIZE
# outputs:
#     bin counts for all possible chrom bins on STDOUT

# get arguments
export BAM_FILE=$1
export CHROM=$2

# get the maximum possible bin number from the chromosome size
# used by count_chrom_bins.pl
CHROM_SIZE=`awk '$1=="'$CHROM'"' ${GENOME_FASTA_SHM}.fai | cut -f2`
export MAX_BIN_I0=$(((${CHROM_SIZE} - 1) / ${BIN_SIZE}))

# if needed, index the bam file
if [ ! -f ${BAM_FILE}.bai ]; then
    samtools index ${BAM_FILE}
fi

# pull all insert endpoints from the bam file
bash $ACTION_DIR/../pull_insert_endpoints.sh | 

# ensure than only unique reads pairs are counted, i.e., de-duplicate reads
bedtools groupby -g 1,2 -c 3 -o distinct |

# count insert endpoints per bin
perl $ACTION_DIR/count_chrom_bins.pl
