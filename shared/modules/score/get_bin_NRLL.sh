# action:
#     analyze the insert size level distribution for each bin individually
# input:
#     $FILENAME_PREFIX as arg1
#     $EMISS_PROBS_FILE as arg2
#     $INSERT_SPANS_DIR
#     $GENOME_BINS_BED
#     $MAPPABILITY_SIZE_LEVELS
#     $BIN_SIZE
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     insert size NRLL for all possible composite bins on STDOUT

# get arguments
export FILENAME_PREFIX=$1
export EMISS_PROBS_FILE=$2

# recover the cache of individual insert spans for this sample
# columns are: chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right
# inserts are pre-filtered for all quality, inclusion, and mappability criteria
zcat ${INSERT_SPANS_DIR}/${FILENAME_PREFIX}.ins_dedup_mpp.bed.gz | 
cut -f 1-4 |

# calculate log likelihoods for each bin for the histone- and protamine-associated states
# output is tab-delimited: LL_histone, LL_protamine, nInserts
perl ${MODULES_DIR}/score/get_bin_likelihoods.pl | 

# calculate NRLL value for each bin
awk '{ print $3 == 0 ? "NA" : ($2 - $1) / $3 }'
