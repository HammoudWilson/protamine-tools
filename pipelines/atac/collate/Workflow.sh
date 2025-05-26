#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/bin/set_genome_bin_vars.sh
source ${MODULES_DIR}/mappability/set_mappability_vars.sh

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh
source ${MODULES_DIR}/utilities/shell/create_temp_dir.sh

# put genome fasta in shared memory
echo "copying genome fasta to shared memory"
export GENOME_FASTA_SHM=$SHM_FILE_PREFIX.genome.fa
cp $GENOME_FASTA $GENOME_FASTA_SHM
cp $GENOME_FASTA.fai $GENOME_FASTA_SHM.fai

# ########################
# export GENOME_FASTA_SHM=$GENOME_FASTA

# establish binned read counts and insert size distributions per sample
runWorkflowStep 1 collate collate.sh

# # calculate and analyze all scores/metrics that do not depend on GC bias normalization
# runWorkflowStep 2 score score.sh

# # profile and establish corrections for in Tn5 site bias
# runWorkflowStep 3 collate bias.sh


