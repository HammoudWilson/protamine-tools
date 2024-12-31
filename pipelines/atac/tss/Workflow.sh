#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/tss/set_tss_vars.sh

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh

# map ATAC fragments around active TSSs at base-level resolution
runWorkflowStep 1 tss tss.sh
