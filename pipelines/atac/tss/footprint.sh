#!/bin/bash

# clean up prior files to avoid problems with merging globs
rm -rf ${TASK_DIR}/inserts_bgz/*.bed.bgz
rm -rf ${TASK_DIR}/inserts_bgz/*.bed.bgz.tbi
rm -rf ${TASK_DIR}/inserts_bgz/*.bed.bgz.rds
rm -rf ${TASK_DIR}/inserts_bgz/*.tally.txt

# do the work
Rscript ${ACTION_DIR}/footprint.R
checkPipe
