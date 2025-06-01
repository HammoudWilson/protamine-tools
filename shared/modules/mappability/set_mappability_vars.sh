# action:
#     set environment variables that define kmer mappability thresholds, a.k.a. insert size levels
#     values are hardcoded based on physical realities of ATAC-seq insert sizes in mammalian genomes
# expects:
#     source ${MODULES_DIR}/genome/set_genome_vars.sh
# usage:
#     source $MODULES_DIR/mappability/set_mappability_vars.sh

# set the profiling parameters
#   20 insert-size levels, with more density (smaller intervals) at smaller insert sizes where things change faster
#   280 is the maximum possible merged insert size with 2x150 paired reads, indexed as single merged reads
#   larger inserts are paired ends, will take mappability as equivalent to 300bp
#   expect mappability to show size thresholds, e.g., map_60=true implies map_80=true, etc.
export MAPPABILITY_SIZE_LEVELS="35 40 45 50 55 60 70 80 90 100 120 140 160 180 200 220 240 260 280 300"
export SIZE_LEVEL_MAXIMA="39 44 49 54 59 69 79 89 99 119 139 159 179 199 219 239 259 279 299 650"

# path to the mappability files
export MAPPABILITY_FILE_PREFIX=${GENOME_GENMAP_DIR}/maps/${GENOME}.mappability
