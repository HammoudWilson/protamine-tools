#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/bin/set_genome_bin_vars.sh

# initialize temporary directories
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh  # for storing the genome fasta
source ${MODULES_DIR}/utilities/shell/create_temp_dir.sh # for storing intermediate files

# calculate and analyze bin scores/metrics
runWorkflowStep 1 score score.sh

# clean up
rm -rf ${SHM_DIR_WRK}
rm -rf ${TMP_DIR_WRK}
