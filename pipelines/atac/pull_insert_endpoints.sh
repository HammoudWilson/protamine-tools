# action:
#     prepare bam alignments for bin counting and insert size assessment
# expects:
#     $ENFORCE_INCLUSIONS as arg1
#     $GENOME_INCLUSIONS_BED, as dictated by $ENFORCE_INCLUSIONS
#     $BAM_FILE
#     $CHROM, can be empty to assess all nuclear chromosomes
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     insert endpoints as BED "chrom\tstart0\tend1" on STDOUT
#       non deduplicated yet
#       chroms names from COMPOSITE_GENOME

# parse arguments
ENFORCE_INCLUSIONS=$1

# if needed, set the command to restrict counting to included regions
COMPRESS_TO_BAM="cat"
INCLUDE_COMMAND="cat"
if [ "$ENFORCE_INCLUSIONS" != "" ]; then
    COMPRESS_TO_BAM="samtools view --bam -"
    INCLUDE_COMMAND="bedtools intersect -a - -b ${GENOME_INCLUSIONS_BED}"
fi

# when pulling all chromosomes from the bam file, restrict to nuclear chromosomes
RESTRICT_CHROMS_COMMAND="cat"
if [ "$CHROM" == "" ]; then
    RESTRICT_CHROMS_COMMAND="perl $ACTION_DIR/../filter_chroms.pl"
fi

# pull fully aligned read (pairs), enforcing minimum mapping quality
#   --exclude-flags 1804
#       read unmapped (0x4)
#       mate unmapped (0x8)
#       not primary alignment (0x100)
#       read fails platform/vendor quality checks (0x200)
#       read is PCR or optical duplicate (0x400)
samtools view --with-header --min-MQ $MIN_MAPQ --exclude-flags 1804 $BAM_FILE $CHROM | 

# apply final FLAG filters, passing only merged reads and proper pairs
# must do this before filtering to inclusions as supplemental alns are required for filtering merged reads
perl $ACTION_DIR/../pull_proper_reads.pl |

# restrict to included regions as requested
$COMPRESS_TO_BAM | 
$INCLUDE_COMMAND |
samtools view - | 

# restrict to nuclear chromosomes as needed
$RESTRICT_CHROMS_COMMAND |

# enforce insert size length filters
# given filtering above, all TLEN will be positive (TLEN is set by pull_proper_reads.pl for merged reads)
# output as BED3 (i.e., chrom, 0-indexed start, 1-indexed end) over entire read pair insert
awk 'BEGIN{OFS="\t"}$9>='$MIN_INSERT_SIZE'&&$9<='$MAX_INSERT_SIZE'{print $3, $4 - 1, $4 + $9 - 1}'
