# action:
#     run ab initio positioned nucleosome scoring on a spermatid stage
# input:
#     $STAGE_TYPE as arg1
#     $STAGE as arg2
# outputs:
#     PENDING tabix-indexed BGZ_FILE in $TASK_DIR/inserts_bgz

# get arguments
export STAGE_TYPE=$1
export STAGE=$2
export BGZ_FILE_PREFIX=${STAGE_TYPE}.${STAGE}
export BGZ_FILE_NAME=${BGZ_FILE_PREFIX}.bed.bgz
export BGZ_FILE=${TASK_DIR}/inserts_bgz/${BGZ_FILE_NAME}

export NFR_FILE=${TASK_DIR}/ab_initio_bgz/${BGZ_FILE_PREFIX}.nfr.bed.bgz

export AB_INITIO=/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/protamine-tools/pipelines/atac/tss/crates/target/debug/ab_initio

# tabix ${BGZ_FILE} chr1:100000000-150000000 | 
# ${AB_INITIO}

zcat ${BGZ_FILE} | 
${AB_INITIO} | 
bgzip -c > ${NFR_FILE}
tabix -p bed ${NFR_FILE}

zcat ${NFR_FILE}
