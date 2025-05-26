#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/tss/set_tss_vars.sh

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh

# map ATAC fragments around active TSSs at base-level resolution
runWorkflowStep 1 tss tss.sh

# further attempt to find positioned nucleosomes distal to active TSSs
runWorkflowStep 2 find_nuc find_nuc.sh

# # find positioned nucleosomes independently of TSSs (but score them based on proximity to TSSs)
# runWorkflowStep 3 ab_initio ab_initio.sh
