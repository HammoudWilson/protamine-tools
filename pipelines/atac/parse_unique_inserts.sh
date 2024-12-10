# action:
#     prepare bam alignments for bin counting and insert size assessment
# expects:
#     $ENFORCE_EXCLUSIONS as arg1
#     $GENOME_EXCLUSIONS_BED and $GENOME_GAPS_FILE, as dictated by $ENFORCE_EXCLUSIONS
#     $BAM_FILE
#     $CHROM, can be empty to assess all chromosomes
#     $BIN_SIZE
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     distinct insert endpoints in format "start1\tend1_1[,end1_2,...]"on STDOUT

# parse arguments
ENFORCE_EXCLUSIONS=$1
EXCLUDE_COMMAND="cat"
if [ "$ENFORCE_EXCLUSIONS" != "" ]; then
    EXCLUDE_COMMAND="bedtools intersect -v -a - -b ${GENOME_EXCLUSIONS_BED}"
fi

# if needed, remove "chr" from $CHROM to match the bam file
RNAME=`samtools view $BAM_FILE | head -n 1 | cut -f 3`
CHROM_WRK=$CHROM
if [ "$CHROM" != "" ]; then
    if [[ $RNAME == "chr"* && $CHROM_WRK != "chr"* ]]; then
        CHROM_WRK="chr"$CHROM_WRK
    elif [[ $RNAME != "chr"* && $CHROM_WRK == "chr"* ]]; then
        CHROM_WRK=${CHROM_WRK:3}
    fi
fi

# pull the forward read only (aligns to the left in the genome) from high-quality, properly-paired reads
#   --require-flags 3
#       read paired (0x1)
#       read mapped in proper pair (0x2)
#   --exclude-flags 3868
#       read unmapped (0x4)
#       mate unmapped (0x8)*
#       read reverse strand (0x10)
#       not primary alignment (0x100)
#       read fails platform/vendor quality checks (0x200)
#       read is PCR or optical duplicate (0x400)
#       supplementary alignment (0x800)
# enforce minimum mapping quality
samtools view --bam --require-flags 3 --exclude-flags 3868 --min-MQ $MIN_MAPQ $BAM_FILE $CHROM_WRK | 

# remove excluded regions as requested
# TODO: this doesn't work propery when CHROM_WRK has been changed, i.e., when GENOME_EXCLUSIONS_BED and BAM_FILE have different chromosome naming conventions
$EXCLUDE_COMMAND |
samtools view |

# enforce insert size length filters
# also, disallow inserts longer than the bin size, so all reads will be counted in at most two bins (one per endpoint)
# given filtering above, all TLEN will be positive
# output is start1,end1 (i.e., 1-indexed start,end) over entire read pair insert
awk 'BEGIN{OFS="\t"}$9>='$MIN_INSERT_SIZE'&&$9<='$MAX_INSERT_SIZE'&&$9<='$BIN_SIZE'{print $4, $4 + $9 - 1}' | 

# ensure than only unique reads pairs are counted, i.e., finish de-duplicating reads
bedtools groupby -g 1 -c 2 -o distinct
