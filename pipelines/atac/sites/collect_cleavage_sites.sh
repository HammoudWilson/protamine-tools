# action:
#     working from a single sample bam file:
#       filter (merged) read pair, i.e., insert, alignments against:
#           minimum mapping quality
#           nuclear chromosomes
#           insert size limits
#           reference exclusions
#       deduplicate the inserts to at most one insert counted per unique genomic span
#       filter deduplicated inserts against genmap mappability by insert size level
#       save the insert spans for atac/collate
#       collect the reference sequence of the genomic span of each allowed insert, with Tn5 flank padding
#       count Tn5 insertion site nonamers
# input:
#     $BAM_FILE as arg1
#     $TMP_FILE_PREFIX as arg2
#     $SAMPLE_NAME as arg3
#     $FILENAME_PREFIX as arg4
#     $GENOME_FASTA_SHM, on /dev/shm for maxmimum lookup performance
#     $GENOME_GAPS_FILE
#     $GENOME_EXCLUSIONS_BED
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $PRIMARY_GENOME
#     $SPIKE_IN_GENOME
#     $MAPPABILITY_FILE_PREFIX
#     $MAPPABILITY_SIZE_LEVELS
#     $TN5_FLANK_PADDING
#     $TN5_KMERS_FILE
#     $TN5_PREFERENCE_POSITIONS
#     $TASK_DIR
# output:
#     insert spans files in $INSERT_SPANS_DIR
#     Tn5 nonamer counts on STDOUT

# get arguments
export BAM_FILE=$1
export TMP_FILE_PREFIX=$2
export SAMPLE_NAME=$3
export FILENAME_PREFIX=$4
export STAGE_COUNT_FILE=${TMP_FILE_PREFIX}.n_ins_stage.txt # number of inserts after each filtering stage, for logging

# pull fully aligned read (pairs), enforcing minimum mapping quality
#   --exclude-flags 1804
#       read unmapped (0x4)
#       mate unmapped (0x8)
#       not primary alignment (0x100)
#       read fails platform/vendor quality checks (0x200)
#       read is PCR or optical duplicate (0x400)
samtools view --min-MQ ${MIN_MAPQ} --exclude-flags 1804 ${BAM_FILE} | 

##############################
# head -n 1000000 | 

# apply final FLAG filters, passing only merged reads and proper pairs
# output has one BED3 line per (merged) read pair span
perl ${ACTION_DIR}/pull_proper_bed3.pl |
tee >(wc -l | awk '{print "'${SAMPLE_NAME}'\tstage_1_proper_pairs\t"$1}' > ${STAGE_COUNT_FILE}) |

# remove any inserts that cross into any excluded region
bedtools intersect -v -a - -b ${GENOME_EXCLUSIONS_BED} <(cut -f2-4 ${GENOME_GAPS_FILE}) | 
tee >(wc -l | awk '{print "'${SAMPLE_NAME}'\tstage_2_included\t"$1}' >> ${STAGE_COUNT_FILE}) |

# deduplicate and analyze insert genome spans
# reject inserts that genmap called unmappable, i.e., that don't start in a mappable run of starts at the insert size level
#   obviously, bowtie2 thought all alignments to be mappable but our bin weights expect genmap mappable inserts
# perl script handles stage count file entry
perl ${ACTION_DIR}/filter_dedup_mappable.pl

# merge the unique filtered inserts at each insert size level
# each row is a unique insert as chrom, start0, end1, insertSizeLevel
# despite this sort merge, the final ins_dedup_mpp.bed.gz insert_spans file is not guaranteed to be sorted (and is not)
#   unclear why, possibly bedtools getfasta does this?
#   at present, no downstream code requires sorted insert spans but caution should be taken
sort --merge -k1,1 -k2,2n -k3,3n ${TMP_FILE_PREFIX}.ins_dedup_mpp_*.txt | 
tee >(wc -l | awk '{print "'${SAMPLE_NAME}'\tstage_4_mappable\t"$1}' >> $STAGE_COUNT_FILE) |

# adjust genomic spans to include the Tn5 cleavage flanks prior to reference sequence extraction
awk -v p=${TN5_FLANK_PADDING} 'BEGIN{OFS="\t"}$2 >= p{
    $2 -= p; 
    $3 += p;
    print $0;
}' | 

# collect the sequence of the entire reference span to which the read aligned
# used downstream to assess cleavage site usage
# output is chrom, start0, end1, insertSizeLevel, ref_seq
bedtools getfasta -bedOut -fi ${GENOME_FASTA_SHM} -bed - |

# count Tn5 cleavage site nonamers
# write a file of all sample insert spans with Tn5 cleavage site nonamers
# write aggregated nonamer usage to STDOUT
perl ${ACTION_DIR}/count_cleavage_sites.pl
