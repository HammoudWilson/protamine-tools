# action:
#     set environment variables for fixed-width genome bins
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     $BIN_SIZE
# usage:
#     source $MODULES_DIR/bin/set_genome_bin_vars.sh

mkdir -p ${GENOME_BINS_DIR}
export GENOME_BINS_PREFIX=${GENOME_BINS_DIR}/${GENOME}
export GENOME_BINS_BED=${GENOME_BINS_PREFIX}.bin_${BIN_SIZE}.bed.gz
