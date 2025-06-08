# action:
#     collect, weight, sort, and index all primary genome inserts for one sample
# input:
#     $STAGE_TYPE as arg1
#     $STAGE as arg2
#     $SAMPLE_NAME as arg3
#     $FILENAME_PREFIX as arg4
#     $BYTES_RAM_PER_SORT as arg5
#     $INSERT_SPANS_DIR
#     $MAPPABILITY_SIZE_LEVELS
#     $PRIMARY_GENOME
#     $TN5_KMERS_FILE
#     $DATA_FILE_PREFIX
#     $GENOME_METADATA_PREFIX
#     $SHM_DIR_WRK
#     $MAX_SITE_WEIGHT
# outputs:
#     tabix-indexed BGZ_FILE in $TASK_DIR/inserts_bgz

# get arguments
export STAGE_TYPE=$1
export STAGE=$2
export SAMPLE_NAME=$3
export FILENAME_PREFIX=$4
export BYTES_RAM_PER_SORT=$5
export BGZ_FILE_PREFIX=${STAGE_TYPE}.${STAGE}.${SAMPLE_NAME}
export BGZ_FILE_NAME=${BGZ_FILE_PREFIX}.bed.bgz
export BGZ_FILE=${TASK_DIR}/inserts_bgz/${BGZ_FILE_NAME}
export TALLY_FILE=${TASK_DIR}/inserts_bgz/${BGZ_FILE_PREFIX}.tally.txt

# recover the cache of individual insert spans for this sample
# columns are: chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right
# inserts are pre-filtered for all quality, inclusion, and mappability criteria
zcat ${INSERT_SPANS_DIR}/${FILENAME_PREFIX}.ins_dedup_mpp.bed.gz | 

# parse inserts into endpoint site weights
# filter to primary genome only and simplify chrom name
perl ${ACTION_DIR}/footprint_sample.pl | 

# sort (not necessarily sorted to this point)
sort -k1,1 -k2,2n --buffer-size=1G --buffer-size ${BYTES_RAM_PER_SORT}"b" |

# bzip and index
bgzip -c > ${BGZ_FILE}
tabix -p bed ${BGZ_FILE}

# return the BGZ_FILE basename and counts
echo -e "bgzFileName\tnInserts\twInserts"
cat ${TALLY_FILE}
