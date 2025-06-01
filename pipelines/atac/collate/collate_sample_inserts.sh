# action:
#     working from the files written by atac/sites and genome/mappability:
#       account for Tn5 cleavage site preference by calculating counting weights
#       count insert sizes vs. GC content over the read span (n_is_gc)
#       report individual inserts with bin and insert size level on STDOUT (ins_bin_isl)
# input:
#     $FILENAME_PREFIX as arg1
#     $TMP_FILE_PREFIX
#     $INSERT_SPANS_DIR
#     $GENOME_FASTA_SHM, on /dev/shm for maxmimum lookup performance
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $BIN_SIZE
#     $MAPPABILITY_SIZE_LEVELS
#     $TN5_KMERS_FILE
#     $GENOME_METADATA_PREFIX
#     $PRIMARY_GENOME
#     $SPIKE_IN_GENOME

# get arguments
export FILENAME_PREFIX=$1

# recover the cache of individual insert spans for this sample
# columns are: chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right
# inserts are pre-filtered for all quality, inclusion, and mappability criteria
zcat ${INSERT_SPANS_DIR}/${FILENAME_PREFIX}.ins_dedup_mpp.bed.gz | 

#################################
# head -n 100000 |

# get the fraction GC of the insert span on reference
bedtools nuc -fi ${GENOME_FASTA_SHM} -bed - |
tail -n +2 | # remove the bedtools nuc header

# extract bedtools nuc GC column and round to 4 decimal places
awk 'BEGIN{OFS="\t"}{ 
    print $1, $2, $3, $4, $5, $6, int($8 * 10000 + 0.5) / 10000 
}' |

# output final inserts with metadata and insert size distributions stratified by GC
perl $ACTION_DIR/collate_sample_inserts.pl
