#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/genome/set_spike_in_vars.sh
source ${MODULES_DIR}/bin/set_genome_bin_vars.sh
source ${MODULES_DIR}/bin/set_spike_in_bin_vars.sh

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh

# establish binned read counts and insert size distributions per sample
runWorkflowStep 1 collate collate.sh

# calculate and analyze all scores/metrics that do not depend on GC bias normalization
runWorkflowStep 2 score score.sh
