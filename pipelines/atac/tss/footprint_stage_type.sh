# action:
#     collect, weight, sort, and index all primary genome inserts for one spermatid stage type
# input:
#     $STAGE_TYPE as arg1
#     $BYTES_RAM_PER_SORT as arg2
#     $SHM_DIR_WRK
# outputs:
#     tabix-indexed BGZ_FILE in $TASK_DIR/inserts_bgz

# get arguments
export STAGE_TYPE=$1
export BYTES_RAM_PER_SORT=$2
export BGZ_FILE_PREFIX=${STAGE_TYPE}
export BGZ_FILE_NAME=${BGZ_FILE_PREFIX}.bed.bgz
export BGZ_FILE=${TASK_DIR}/inserts_bgz/${BGZ_FILE_NAME}
export TALLY_GLOB=${TASK_DIR}/inserts_bgz/${BGZ_FILE_PREFIX}.*.*.tally.txt
export INPUT_GLOB=${TASK_DIR}/inserts_bgz/${BGZ_FILE_PREFIX}.*.*.bed.bgz # does not match stage-level files, only sample-level files

# sort input files together
zcat ${INPUT_GLOB} | 
sort -k1,1 -k2,2n --buffer-size=${BYTES_RAM_PER_SORT}b | 

# bzip and index
bgzip -c > ${BGZ_FILE}
tabix -p bed ${BGZ_FILE}

# return the BGZ_FILE basename and counts
echo -e "bgzFileName\tnInserts\twInserts"
cat ${TALLY_GLOB} |
awk -v bgz=${BGZ_FILE_NAME} 'BEGIN {OFS="\t"} {print bgz, $2, $3}' | 
bedtools groupby -g 1 -c 2,3 -o sum,sum
