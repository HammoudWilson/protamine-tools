#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/genome/set_spike_in_vars.sh
source ${MODULES_DIR}/bin/set_genome_bin_vars.sh
source ${MODULES_DIR}/bin/set_spike_in_bin_vars.sh

# establish binned read counts and insert size distributions per sample
runWorkflowStep 1 collate collate.sh
