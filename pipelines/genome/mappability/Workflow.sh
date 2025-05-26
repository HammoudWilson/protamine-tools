#!/bin/bash

# set the output directories
source $MODULES_DIR/genome/set_genome_vars.sh
source $MODULES_DIR/mappability/set_mappability_vars.sh
mkdir -p $GENOME_GENMAP_DIR

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh

# put genome fasta in shared memory
echo "copying genome fasta to shared memory"
export GENOME_FASTA_SHM=$SHM_FILE_PREFIX.genome.fa
cp $GENOME_FASTA $GENOME_FASTA_SHM
cp $GENOME_FASTA.fai $GENOME_FASTA_SHM.fai

# create a genome bins file at the highest resolution with exclusion and GC metadata
runWorkflowStep 1 mappability mappability.sh
