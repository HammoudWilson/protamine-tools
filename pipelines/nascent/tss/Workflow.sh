#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/tss/set_tss_vars.sh

# create a genome bins file at the highest resolution with exclusion and GC metadata
runWorkflowStep 1 tss tss.sh
