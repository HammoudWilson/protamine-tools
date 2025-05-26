# action:
#     set environment variables that define
#       allowed insert size ranges
#       kmer mappability thresholds
#     values are hardcoded based on physical realities of ATAC-seq insert sizes mammalian genomes
# expects:
#     source ${MODULES_DIR}/genome/set_genome_vars.sh
# usage:
#     source $MODULES_DIR/mappability/set_mappability_vars.sh

# set the allowed insert size ranges
export MIN_INSERT_SIZE=35   # smaller inserts are poorly mappable and show increasing GC bias
export MAX_INSERT_SIZE=650  # ~corresponds to the upper limit of tri-nucleosome size fragments

# set the profiling parameters
#   20 insert-size levels, with more density at smaller inserts
#   280 is the maximum merged insert size with 2x150 paired reads, indexed as single merged reads
#   larger inserts are paired ends, will take mappability as equivalent to 300bp
#   expect mappability to show size thresholds, e.g., map_60=true implies map_80=true, etc.
export MAPPABILITY_KMER_LENGTHS="35 40 45 50 55 60 70 80 90 100 120 140 160 180 200 220 240 260 280 300" 

# path to the mappability files
export MAPPABILITY_FILE_PREFIX=${GENOME_GENMAP_DIR}/maps/${GENOME}.mappability

# the positions of Tn5 cleavage sites that most impact insert recovery
# empirical logo plots reveal the following positions as having greatest Tn5 preference
#   first  index is with respect to the padded insert span
#   second index is with respect to the central position of the 9-mer overhang
#  0 -6    +++    # strength of preference -/+/++/+++/++++
#  1 -5    ++
# -----------     # cleaved bond on top strand
#  2 -4    ++++   # most informative position; bases upstream of insert 5' end have more impact
#  3 -3    + 
#  4 -2    -      # little to no preference
#  5 -1    +++
#  6  0    ++     # central position of the 9-mer overhang
#  7  1    +++
#  8  2    -
#  9  3    -
# 10  4    +++
# -----------     # cleaved bond on bottom strand
# 11  5    +
# 12  6    ++
# note that the Tn5 preference is not symmetric
# the map above is for the top strand, or the reverse complement of the bottom strand
#    0  v  (carats denote insert span ends after 5' overhang filling)
# ---32|410323003 12---
# ---21 300323014|23---
#               ^  0
# for efficiency, restrict site matching to the three trinucleotide cluster
# yields a total of 9 preferred positions for 262,144 kmers
#    0  v
# ---**|*--***--* **---
# ---** *--***--*|**---
#               ^  0
export TN5_FLANK_PADDING=2
export TN5_PREFERENCE_POSITIONS="0 1 2 5 6 7 10 11 12"
export N_TN5_PREFERENCE_POSITIONS=9
export TN5_PREFERENCE_SPAN=13
export TN5_KMERS_FILE=${TASK_DIR}/Tn5_kmers.txt
if [ ! -f $TN5_KMERS_FILE ]; then
    echo "writing file of all possible Tn5 $N_TN5_PREFERENCE_POSITIONS-kmers"
    perl ${MODULES_DIR}/mappability/create_Tn5_kmers.pl
fi
