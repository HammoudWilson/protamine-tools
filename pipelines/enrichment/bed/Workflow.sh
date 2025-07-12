#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_composite_vars.sh
source ${MODULES_DIR}/genome/set_genome_vars.sh
export REGIONS_BED_DIR=${TASK_DIR}/regions_bed
mkdir -p ${REGIONS_BED_DIR}

# map regions BED files from the list in regions_types.txt
# only types in that file will be created, regardless the scripts found in the subfolders
runWorkflowStep 1 make_bed make_bed.sh
