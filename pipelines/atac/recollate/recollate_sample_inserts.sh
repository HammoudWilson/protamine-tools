# action:
#     working from the files written by atac/sites and genome/mappability:
#       count insert sizes vs. GC content over the read span (n_is_gc)
#       report individual inserts with bin and insert size level on STDOUT (ins_bin_isl)
#       working with observed counts and primary genome only
# input:
#     $FILENAME_PREFIX as arg1
#     $PEAK_REGIONS_BED as arg2
#     $TMP_FILE_PREFIX
#     $INSERT_SPANS_DIR
#     $GENOME_FASTA_SHM, on /dev/shm for maxmimum lookup performance
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $BIN_SIZE
#     $MAPPABILITY_SIZE_LEVELS
#     $GENOME_METADATA_PREFIX
#     $PRIMARY_GENOME

# get arguments
export FILENAME_PREFIX=$1
export PEAK_REGIONS_BED=$2

# recover the cache of individual insert spans for this sample
# columns are: chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right
# inserts are pre-filtered for all quality, inclusion, and mappability criteria
zcat ${INSERT_SPANS_DIR}/${FILENAME_PREFIX}.ins_dedup_mpp.bed.gz | 

# get the fraction GC of the insert span on reference
bedtools nuc -fi ${GENOME_FASTA_SHM} -bed - |
tail -n +2 | # remove the bedtools nuc header

# extract bedtools nuc GC column and round to 4 decimal places
# drop unused Tn5 information
awk 'BEGIN{OFS="\t"}{ 
    print $1, $2, $3, $4, int($8 * 10000 + 0.5) / 10000 
}' |

# flag inserts that overlap a peak region
bedtools intersect -c -a - -b ${PEAK_REGIONS_BED} |

# output final inserts with metadata and insert size distributions stratified by GC
perl $ACTION_DIR/recollate_sample_inserts.pl
