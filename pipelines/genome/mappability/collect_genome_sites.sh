# action:
#     count all possible discontiguous Tn5 nonamers within genome mappable regions at given insert size level
# input:
#     $MPP_FILE as arg1
#     $INSERT_SIZE_LEVEL as arg2
#     $MAX_INSERT_SIZE as arg3
#     $GENOME_FASTA_SHM, on /dev/shm for maxmimum lookup performance
#     $PRIMARY_GENOME
#     $SPIKE_IN_GENOME
#     $TN5_FLANK_PADDING
#     $TN5_KMERS_FILE
#     $TN5_PREFERENCE_POSITIONS
#     $TMP_FILE_PREFIX
# output:
#     Tn5 nonamer counts on STDOUT

# get arguments
export MPP_FILE=$1
export INSERT_SIZE_LEVEL=$2
export MAX_INSERT_SIZE=$3
MASKED_FASTA=${TMP_FILE_PREFIX}.k_${INSERT_SIZE_LEVEL}.masked.fa

# collect all runs of mappable start positions at the insert size level
zcat ${MPP_FILE} |
cut -f 1-3 |

# pad the run spans to match all bases that might contribute to Tn5 cleavage sites in kept inserts
awk -v p=${TN5_FLANK_PADDING} -v mis=${MAX_INSERT_SIZE} 'BEGIN {OFS="\t"}$2 >= p{
    $2 -= p;
    $3 += (mis - 1 + p);
    print $0;
}' |

# remove mappable regions on non-target chromosomes
# those chroms get complemented to fully masked regions
# count_genome_sites.pl bypasses them in any case, but this line suppresses bedtools warnings
grep -v -e "_" -e "chrM" -e "chrEBV" | 

# get the complement of the bases above, i.e., bases that could never contribute to mappable inserts
# implicitly merges overlapping padded mappable regions, which will be scanned end to end for Tn5 cleavage sites
# turns out, there are still rare occasions when code above exceeds the end of a chromosome, so suppress warnings
bedtools complement -i - -g ${GENOME_FASTA_SHM}.fai 2>/dev/null |

# mask the complement regions as N bases to prevent Tn5 cleavage site counting in unmappable bases
# I believe the masked genome must be written to file first; this step is fast enough regardless
bedtools maskfasta -bed - -fi ${GENOME_FASTA_SHM} -fo ${MASKED_FASTA}

# count Tn5 cleavage site nonamers using the masked fasta
cat ${MASKED_FASTA} |
perl ${ACTION_DIR}/count_genome_sites.pl
