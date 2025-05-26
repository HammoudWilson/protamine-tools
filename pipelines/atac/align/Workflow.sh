#!/bin/bash

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_temp_dir.sh

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source $MODULES_DIR/mappability/set_mappability_vars.sh

# create the bowtie2 index if needed
runWorkflowStep 1 index index.sh

# use fastp to trim, filter and merge the reads
# use Bowtie2 to align the reads to the genome, one sample at a time
runWorkflowStep 2 align align.sh

# clean up 
rm -rf $TMP_DIR_WRK
