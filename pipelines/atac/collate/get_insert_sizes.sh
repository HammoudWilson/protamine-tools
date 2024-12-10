# action:
#     calculate insert size distributions in a single bam file
# input:
#     $BAM_FILE as arg1
#     $BIN_SIZE, for excluding overly large inserts to match bin counts
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     insert size distribution on STDOUT

# get arguments
export BAM_FILE=$1
export CHROM="" # thus, this script acts over all chromosomes

# pull all distinct insert endpoints from the bam file
# since final output is aggregated, enforce exclusions here rather than downstream
bash $ACTION_DIR/../parse_unique_inserts.sh ENFORCE_EXCLUSIONS | 

# count insert sizes
perl $ACTION_DIR/get_insert_sizes.pl
