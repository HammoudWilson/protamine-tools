# action:
#     collect read, i.e., insert, endpoints in TSS regions for a single sample
# input:
#     $FILENAME_PREFIX as arg1
#     $TSS_FILE as arg2
#     $INSERT_SPANS_DIR
#     $GENOME_FASTA
# outputs:
#     one line per insert in a TSS region on STDOUT (see Perl)

# get arguments
export FILENAME_PREFIX=$1
export TSS_FILE=$2

# recover the cache of individual insert spans for this sample
# columns are: chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right
# inserts are pre-filtered for all quality, inclusion, and mappability criteria
zcat ${INSERT_SPANS_DIR}/${FILENAME_PREFIX}.ins_dedup_mpp.bed.gz | 
cut -f 1-3 |

# filter de-duplicated inserts against active TSS regions
# -wo only keeps inserts that overlap TSS regions, and reports each TSS overlap
# i.e., inserts may appear for each of two closely spaced TSSs
# NB: this insersect is not stranded; ATAC-seq is not strand-specific
# TSS strandedness is accounted for during distance calculation below
# sort order of ins_dedup_mpp.bed.gz is not guaranteed
bedtools intersect -wo -a - -b <(
    zcat ${TSS_FILE} | 
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,NR}' # add a unique ID to each TSS
) |

# calculate TSS distances for eventual plotting
perl ${ACTION_DIR}/get_tss_inserts.pl
