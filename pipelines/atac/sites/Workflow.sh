#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/insert/set_insert_vars.sh
source ${MODULES_DIR}/mappability/set_mappability_vars.sh
source ${MODULES_DIR}/tn5/set_tn5_site_vars.sh
mkdir -p ${INSERT_SPANS_DIR}

# initialize temporary directories
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh  # for storing the genome fasta
source ${MODULES_DIR}/utilities/shell/create_temp_dir.sh # for storing intermediate insert size files

# copy genome fasta into shared memory for speed
echo "copying genome fasta to shared memory"
export GENOME_FASTA_SHM=${SHM_FILE_PREFIX}.genome.fa
cp ${GENOME_FASTA} ${GENOME_FASTA_SHM}
cp ${GENOME_FASTA}.fai ${GENOME_FASTA_SHM}.fai

# establish binned read counts and insert size distributions per sample
runWorkflowStep 1 sites sites.sh
