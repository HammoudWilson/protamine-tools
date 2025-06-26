# action:
#     set environment variables for creating a composite genome
# expects:
#     genome/download has been completed on both the primary and spike-in genomes
#     all manual exclusion files have also been downloaded to create both verions of:
#         GENOME_FASTA
#         GENOME_GAPS_FILE
#         GENOME_EXCLUSIONS_BED
# usage:
#     source $MODULES_DIR/genome/set_composite_vars.sh

# save GENOME and GENOME_DIR variables
export GENOME_TMP=${GENOME}
export GENOME_DIR_TMP=${GENOME_DIR}

# use set_genome_vars to parse the primary genome
export GENOME=${PRIMARY_GENOME}
export GENOME_DIR=${PRIMARY_GENOME_DIR}
source ${MODULES_DIR}/genome/set_genome_vars.sh
export PRIMARY_GENOME_FASTA=${GENOME_FASTA}
export PRIMARY_GENOME_GAPS_FILE=${GENOME_GAPS_FILE}
export PRIMARY_GENOME_EXCLUSIONS_BED=${GENOME_EXCLUSIONS_BED}
export PRIMARY_GENOME_GENES_BED=${GENES_BED}

# use set_genome_vars to parse the spike-in genome (genes bed not used)
export GENOME=${SPIKE_IN_GENOME}
export GENOME_DIR=${SPIKE_IN_GENOME_DIR}
source ${MODULES_DIR}/genome/set_genome_vars.sh
export SPIKE_IN_GENOME_FASTA=${GENOME_FASTA}
export SPIKE_IN_GENOME_GAPS_FILE=${GENOME_GAPS_FILE}
export SPIKE_IN_GENOME_EXCLUSIONS_BED=${GENOME_EXCLUSIONS_BED}

# restore the original genome variables
export GENOME=${GENOME_TMP}
export GENOME_DIR=${GENOME_DIR_TMP}
