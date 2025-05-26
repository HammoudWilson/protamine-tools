# action:
#     working from a single sample bam file:
#       filter (merged) read pair, i.e., insert, alignments against:
#           minimum mapping quality
#           nuclear chromosomes
#           insert size limits
#           reference exclusions
#       deduplicate the inserts to at most one insert counted per unique genomic span
#       filter deduplicated inserts against genmap mappability by insert size level
#       collect the reference sequence of the genomic span of each allowed insert, with Tn5 flank padding
#       count insert sizes vs. GC content over the read span (n_is_gc)
#       count Tn5 insertion site kmers (n_tn5_smp)
#       report individual inserts with bin and insert size level on STDOUT (ins_bin_isl)
# input:
#     $BAM_FILE as arg1
#     $TMP_FILE_PREFIX as arg2
#     $SAMPLE_LABEL as arg3
#     $GENOME_FASTA_SHM, on /dev/shm for maxmimum lookup performance
#     $GENOME_GAPS_FILE
#     $GENOME_EXCLUSIONS_BED
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $PRIMARY_GENOME
#     $SPIKE_IN_GENOME
#     $TN5_FLANK_PADDING
#     $TN5_KMERS_FILE
#     $TN5_PREFERENCE_POSITIONS

# get arguments
export BAM_FILE=$1
export TMP_FILE_PREFIX=$2
export SAMPLE_LABEL=$3
export STAGE_COUNT_FILE=${TMP_FILE_PREFIX}.n_ins_stage.txt # number of inserts after each filtering stage

# pull fully aligned read (pairs), enforcing minimum mapping quality
#   --exclude-flags 1804
#       read unmapped (0x4)
#       mate unmapped (0x8)
#       not primary alignment (0x100)
#       read fails platform/vendor quality checks (0x200)
#       read is PCR or optical duplicate (0x400)
samtools view --min-MQ $MIN_MAPQ --exclude-flags 1804 $BAM_FILE | 

################################
# head -n 1000000 | 

# apply final FLAG filters, passing only merged reads and proper pairs
# output has one BED3 line per (merged) read pair span
perl $ACTION_DIR/pull_proper_bed3.pl |
tee >(wc -l | awk '{print "'${SAMPLE_LABEL}'\tstage_1_proper_pairs\t"$1}' > $STAGE_COUNT_FILE) |

# remove any inserts that cross into any excluded region
bedtools intersect -v -a - -b ${GENOME_EXCLUSIONS_BED} <(cut -f2-4 ${GENOME_GAPS_FILE}) | 
tee >(wc -l | awk '{print "'${SAMPLE_LABEL}'\tstage_2_included\t"$1}' >> $STAGE_COUNT_FILE) |

# deduplicate and analyze insert genome spans
# reject any inserts, i.e., kmers, that genmap called unmappable
#   obviously, bowtie2 thought all alignemnts to be mappable
#   but our bin weights only include genmap mappable kmers
# perl script handle stage count file entry
perl $ACTION_DIR/collate_dedup_mappable.pl

# merge the unique filtered inserts at each insert size level
# each row is a unique insert as chrom, start0, end1, insertSizeLevel
sort --merge -k1,1 -k2,2n -k3,3n $TMP_FILE_PREFIX.ins_dedup_mpp_*.txt | 
tee >(wc -l | awk '{print "'${SAMPLE_LABEL}'\tstage_4_mappable\t"$1}' >> $STAGE_COUNT_FILE) |

# adjust genomic spans to include the Tn5 cleavage flanks prior to reference sequence extraction
awk 'BEGIN{OFS="\t"}$2 >= '$TN5_FLANK_PADDING'{
    print $1, $2 - '$TN5_FLANK_PADDING', $3 + '$TN5_FLANK_PADDING', $4
}' | 

# collect the sequence of the entire reference span to which the read aligned
# used downstream to assess cleavage kmers and GC content
# output is BED3, insert_size_level, ref_seq
bedtools getfasta -bedOut -fi ${GENOME_FASTA_SHM} -bed - |

# output final inserts with metadata and various temporary files, including:
#    insert size distributions stratified by GC
#    observed Tn5 kmer preferences
perl $ACTION_DIR/collate_count_inserts.pl
