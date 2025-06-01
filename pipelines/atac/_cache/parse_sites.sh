# action:
#     prepare bam alignments on all nuclear chroms for insert site assessment
# expects:
#     $ENFORCE_EXCLUSIONS as arg1
#     $EXCLUSIONS_BED, as dictated by $ENFORCE_EXCLUSIONS
#     $BAM_FILE
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $FA_FILE
# outputs:
#     19-mer site sequences on STDOUT
#     these are not de-duplicated, but the impact of duplicates should be minimal and average out

# parse arguments
ENFORCE_EXCLUSIONS=$1

# if needed, set the command to exclude regions from counting
EXCLUDE_COMMAND="cat"
EXCLUSIONS_BED_TMP=$(mktemp)
if [ "$ENFORCE_EXCLUSIONS" != "" ]; then
    EXCLUDE_COMMAND="bedtools intersect -v -a - -b ${EXCLUSIONS_BED_TMP}"

    # if needed, remove "chr" from $EXCLUSIONS_BED to match the bam file
    RNAME=`samtools view $BAM_FILE | head -n 1 | cut -f 3`
    if [[ $RNAME == "chr"* ]]; then
        cp ${EXCLUSIONS_BED} ${EXCLUSIONS_BED_TMP}
    else
        cat ${EXCLUSIONS_BED} | sed 's/^chr//' > ${EXCLUSIONS_BED_TMP}
    fi
fi

# when pulling all chromosome from the bam file, restrict to nuclear chromosomes
RESTRICT_CHROMS_COMMAND="perl $ACTION_DIR/../filter_chroms.pl"

# pull the forward read only (aligns to the left in the genome) from high-quality, properly-paired reads
#   --require-flags 3
#       read paired (0x1)
#       read mapped in proper pair (0x2)
#   --exclude-flags 3868
#       read unmapped (0x4)
#       mate unmapped (0x8)*
#       read reverse strand (0x10)  # for now at least, continue to pull the 5'-most base of forward reads
#       not primary alignment (0x100)
#       read fails platform/vendor quality checks (0x200)
#       read is PCR or optical duplicate (0x400)
#       supplementary alignment (0x800)
# enforce minimum mapping quality

# support mutiple input BAM files, thus potentially concatenating multiple files on the same chromosome
# caller must must be able to account for the fact the alignments are no longer sorted if there are multiple input files
for BAM_FILE_WRK in $BAM_FILE; do 

    samtools view --bam --require-flags 3 --exclude-flags 3868 --min-MQ $MIN_MAPQ $BAM_FILE_WRK | 

    # remove excluded regions as requested
    $EXCLUDE_COMMAND |

    # restrict to nuclear chromosomes as needed
    samtools view |
    $RESTRICT_CHROMS_COMMAND |

    # enforce insert size length filters
    # given filtering above, all TLEN will be positive
    # exclude clipped 5' ends
    # output is BED3, chrom,start0,end1
    awk 'BEGIN{OFS="\t"}$9>='$MIN_INSERT_SIZE'&&$9<='$MAX_INSERT_SIZE'&&$6~/^[0-9]+M/{print $3, $4 - 5 - 1, $4 + 13}' | 
    awk '$2>=0' | 

    # pull the 19-mer sequence from the reference genome
    bedtools getfasta -bedOut -fi $FA_FILE -bed - | 

    # parse to uppercase ACGT only
    cut -f 4 | 
    tr '[:lower:]' '[:upper:]' | 
    grep -E '^[ACGT]{19}$'

done

# clean up tmp file
rm -f $EXCLUSIONS_BED_TMP
