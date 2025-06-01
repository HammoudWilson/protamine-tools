# action:
#     set environment variables that define Tn5 cleavage site discontiguous nonamers
#     the structure is derived from empirical observation of Tn5 cleavage site preferences
# expects:
#     nothing
# usage:
#     source $MODULES_DIR/tn5/set_tn5_site_vars.sh

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

# note that the Tn5 preference is not fully symmetric or palindromic
# the map above is for the top strand, or the reverse complement of the bottom strand
#    0  v  (carats denote insert span ends after 5' overhang filling)
# ---32|410323003 12---
# ---21 300323014|23---
#               ^  0

# for efficiency, restrict site matching to the three trinucleotide clusters
# yields a total of 9 preferred positions for 262,144 kmers
#    0  v
# ---**|*--***--* **---  where * is a base used in a Tn5 cleavage site nonamer
# ---** *--***--*|**---
#               ^  0

# compressed nonamer as recorded in TN5_KMERS_FILE
# *********
# where bases are always recorded from the top strand of the reference genome

# set variable that define discontiguous Tn5 cleavage site nonamers
export TN5_FLANK_PADDING=2
export TN5_PREFERENCE_POSITIONS="0 1 2 5 6 7 10 11 12"
export N_TN5_PREFERENCE_POSITIONS=9
export TN5_PREFERENCE_SPAN=13

# check and create the list of all possible nonamers (frequencies tabulated later)
export TN5_KMERS_FILE=${TASK_DIR}/tn5_kmers.txt
if [ ! -f $TN5_KMERS_FILE ]; then
    echo "writing file of all possible Tn5 $N_TN5_PREFERENCE_POSITIONS-kmers"
    perl ${MODULES_DIR}/tn5/create_Tn5_kmers.pl
fi

