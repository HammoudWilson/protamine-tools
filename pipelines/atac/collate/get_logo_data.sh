# action:
#     create a file with a collection of 19-mer insert size sequences
#     calculate base utilization at each 19-mer position
# input:
#     $BAM_FILE as arg1
#     $GENOME as arg2
#     $FA_FILE as arg3
#     $EXCLUSIONS_BED as arg4
#     $FILE_PREFIX as arg5
#     $STAGING_ORDER as arg6
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
# outputs:
#     insert size distribution on STDOUT

# get arguments
export BAM_FILE=$1
export GENOME=$2
export FA_FILE=$3
export EXCLUSIONS_BED=$4
export SAMPLE_NAME=$5
export STAGING_ORDER=$6
export CHROM="" # thus, this script acts over all chromosomes
export SAMPLE_SITES_FILE=$SAMPLE_SITES_DIR/$SAMPLE_NAME.$GENOME.sites.bed.gz

# pull a representative set of sites from the BAM file
bash $ACTION_DIR/../parse_sites.sh ENFORCE_EXCLUSIONS | 

# head -n 1000000 | 

# extract logo data
perl $ACTION_DIR/extract_logo_data.pl
