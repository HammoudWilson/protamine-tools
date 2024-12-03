#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/bin/set_bin_vars.sh

# create binned read counts per sample
runWorkflowStep 1 bin bin.sh
