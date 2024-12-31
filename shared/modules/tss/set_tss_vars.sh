# action:
#     set environment variables for tss file access
# expects:
#     GENOME
#     GENCODE_RELEASE
# usage:
#     source $MODULES_DIR/tss/set_tss_vars.sh

# active TSS file
export ACTIVE_TSS_BED=${DATA_FILE_PREFIX}.${GENOME}.gencode.${GENCODE_RELEASE}.active.tss.bed.gz
