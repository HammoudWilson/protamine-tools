#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_composite_vars.sh
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/insert/set_insert_vars.sh
source ${MODULES_DIR}/mappability/set_mappability_vars.sh
source ${MODULES_DIR}/tn5/set_tn5_site_vars.sh
source ${MODULES_DIR}/tss/set_tss_vars.sh
mkdir -p ${TASK_DIR}/inserts_bgz
# mkdir -p ${TASK_DIR}/ab_initio_bgz

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh

# map ATAC fragments around active TSSs at base-level resolution
runWorkflowStep 1 tss tss.sh

# further attempt to find positioned nucleosomes distal to active TSSs
runWorkflowStep 2 find_nuc find_nuc.sh

# create tabix-indexed files over primary genome inserts for footprint visualization
runWorkflowStep 3 footprint footprint.sh

# find positioned nucleosomes independently of TSSs (but score them based on proximity to TSSs)
runWorkflowStep 4 ab_initio ab_initio.sh
