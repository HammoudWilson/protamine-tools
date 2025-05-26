# action:
#     calculate nuc (nucleosome fragment coverage) and ends (endpoint count) for all bins on a chromosome
# input:
#     $BAM_FILE     as arg1
#     $CHROM        as arg2
#     $AI_BIN_SIZE  as arg3
#     $MIN_NUC_SIZE as arg4
#     $MAX_NUC_SIZE as arg5
#     $FAI_FILE     as arg6
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     one line per chromosome bin with columns nuc,ends

# get arguments
export BAM_FILE=$1 # can be multiple BAM files
export CHROM=$2
export AI_BIN_SIZE=$3
export MIN_NUC_SIZE=$4
export MAX_NUC_SIZE=$5
export FAI_FILE=$6
export BIN_SIZE=$MAX_INSERT_SIZE # thus, BIN_SIZE does nothing additional in pull_insert_endpoints.sh
export ENFORCE_EXCLUSIONS=1

# get the maximum possible bin number from the chromosome size
# used by score_ab_initio.pl
CHROM_SIZE=`awk '$1=="'$CHROM'"' $FAI_FILE | cut -f2`
export MAX_BIN_I0=$((($CHROM_SIZE - 1) / $AI_BIN_SIZE))

# pull all distinct insert endpoints from the bam file
bash $ACTION_DIR/../pull_insert_endpoints.sh | 

# unpack to sorted BED3
perl $ACTION_DIR/unpack_chrom_inserts.pl | 

# calculate TSS distances for eventual plotting
perl $ACTION_DIR/score_ab_initio.pl
