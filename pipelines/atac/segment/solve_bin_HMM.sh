# action:
#     collect the insert sizes for each bin individually
# input:
#     $BAM_FILE         as arg1
#     $CHROM            as arg2
#     $FAI_FILE         as arg3
#     $EMISS_PROBS_FILE as arg4
#     $PERSISTENCE      as arg5
#     $BIN_SIZE
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     concatenated insert sizes for all possible chrom bins on STDOUT

# get arguments
export BAM_FILE=$1
export CHROM=$2
export FAI_FILE=$3
export EMISS_PROBS_FILE=$4
export PERSISTENCE=$5

# get the maximum possible bin number from the chromosome size
# used by get_bin_likelihoods.pl
CHROM_SIZE=`awk '$1=="'$CHROM'"' $FAI_FILE | cut -f2`
export MAX_BIN_I0=$((($CHROM_SIZE - 1) / $BIN_SIZE))

# pull all distinct insert endpoints from the bam file
bash $ACTION_DIR/../parse_unique_inserts.sh | 

# concatenate unique insert endpoints per bin
perl $ACTION_DIR/get_bin_likelihoods.pl | 

# run the HMM per chromosome
perl $ACTION_DIR/segment.pl
