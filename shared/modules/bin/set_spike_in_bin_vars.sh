# action:
#     set environment variables for fixed-width spike-in-genome bins
# expects:
#     source $MODULES_DIR/genome/set_spike_in_vars.sh
#     $BIN_SIZE
# usage:
#     source $MODULES_DIR/bin/set_spike_in_bin_vars.sh

mkdir -p ${SPIKE_IN_BINS_DIR}
export SPIKE_IN_BINS_PREFIX=${SPIKE_IN_BINS_DIR}/${SPIKE_IN_GENOME}
export SPIKE_IN_BINS_BED=${SPIKE_IN_BINS_PREFIX}.bin_${BIN_SIZE}.bed.gz
