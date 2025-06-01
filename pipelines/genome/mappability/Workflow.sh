#!/bin/bash

# set the output directories
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/mappability/set_mappability_vars.sh
source ${MODULES_DIR}/tn5/set_tn5_site_vars.sh
mkdir -p ${GENOME_GENMAP_DIR}

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh
source ${MODULES_DIR}/utilities/shell/create_temp_dir.sh # for storing intermediate masked genomes

# copy genome fasta into shared memory for speed
echo "copying genome fasta to shared memory"
export GENOME_FASTA_SHM=${SHM_FILE_PREFIX}.genome.fa
cp ${GENOME_FASTA} ${GENOME_FASTA_SHM}
cp ${GENOME_FASTA}.fai ${GENOME_FASTA_SHM}.fai

# profile the composite genome for mappability at different size levels with GC metadata
runWorkflowStep 1 mappability mappability.sh

# count Tn5 site nonamers in the mappable genome regions by insert size level
runWorkflowStep 2 tn5_sites sites.sh

# clean up
rm -rf ${SHM_DIR_WRK}
rm -rf ${TMP_DIR_WRK}
