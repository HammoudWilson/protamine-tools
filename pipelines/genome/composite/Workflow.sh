#!/bin/bash
export CREATING_COMPOSITE=true

# initialize shared variables
source ${MODULES_DIR}/genome/set_composite_vars.sh
source ${MODULES_DIR}/genome/set_genome_vars.sh # must come after set_composite_vars.sh

# set up the output directories
mkdir -p ${GENOME_DIR}/aligner-indices
mkdir -p ${GENOME_DIR}/annotations
mkdir -p ${GENOME_DIR}/bins
mkdir -p ${GENOME_DIR}/metadata

# create a composite genome FASTQ and metadata files
runWorkflowStep 1 composite composite.sh
