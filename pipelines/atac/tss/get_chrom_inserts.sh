# action:
#     collect read, i.e., insert, endpoints in TSS regions for a single chrom in a single bam file
# input:
#     $BAM_FILE     as arg1
#     $CHROM        as arg2
#     $TSS_FILE     as arg3
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     one line per insert in a TSS region on STDOUT (see Perl)

# get arguments
export BAM_FILE=$1
export CHROM=$2
export TSS_FILE=$3
export BIN_SIZE=$MAX_INSERT_SIZE # as required by parse_unique_inserts.sh

# pull all distinct insert endpoints from the bam file
bash $ACTION_DIR/../parse_unique_inserts.sh | 

# unpack to sorted BED3
perl $ACTION_DIR/unpack_chrom_inserts.pl | 

# filter de-duplicated inserts against active TSS regions
# -wo only keeps inserts that overlap TSS regions, and reports each TSS overlap
# i.e., inserts may appear for each of two closely spaced TSSs
# NB: this insersect is not stranded; ATAC-seq is not strand-specific
# TSS strandedness is accounted for during distance calculation below
bedtools intersect -wo -sorted -a - -b <(
    zcat $TSS_FILE | 
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,NR}' | # add a unique ID to each TSS
    awk '$1=="'$CHROM'"'
) |

# calculate TSS distances for eventual plotting
perl $ACTION_DIR/get_chrom_inserts.pl
