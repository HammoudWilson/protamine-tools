# action:
#     set environment variables that define allowed insert size ranges
#     values are hardcoded based on physical realities of ATAC-seq insert sizes in mammalian genomes
# expects:
#     nothing
# usage:
#     source $MODULES_DIR/insert/set_insert_vars.sh

# set the allowed insert size ranges
export MIN_INSERT_SIZE=35   # smaller inserts are poorly mappable and show increasing GC bias
export MAX_INSERT_SIZE=650  # corresponds to the ~upper limit of tri-nucleosome size fragments

# set the path to the insert span files
export INSERT_SPANS_DIR=${TASK_DIR}/insert_spans
