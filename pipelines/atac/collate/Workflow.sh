#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/bin/set_genome_bin_vars.sh
source ${MODULES_DIR}/insert/set_insert_vars.sh
source ${MODULES_DIR}/mappability/set_mappability_vars.sh
source ${MODULES_DIR}/tn5/set_tn5_site_vars.sh

# initialize temporary directories
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh  # for storing the genome fasta
source ${MODULES_DIR}/utilities/shell/create_temp_dir.sh # for storing intermediate files

# copy genome fasta into shared memory for speed
if [ "${PUSHING_TO_SERVER}" = "" ]; then
    echo "copying genome fasta to shared memory"
    export GENOME_FASTA_SHM=${SHM_FILE_PREFIX}.genome.fa
    cp ${GENOME_FASTA} ${GENOME_FASTA_SHM}
    cp ${GENOME_FASTA}.fai ${GENOME_FASTA_SHM}.fai
fi

# establish binned read counts and insert size distributions per sample
runWorkflowStep 1 collate collate.sh

# calculate and analyze bin scores/metrics
runWorkflowStep 2 score ${MODULES_DIR}/score/score.sh

# DEV USE
# # profile and establish corrections for in Tn5 site bias
# runWorkflowStep 3 collate bias.sh

# clean up
rm -rf ${SHM_DIR_WRK}
rm -rf ${TMP_DIR_WRK}
