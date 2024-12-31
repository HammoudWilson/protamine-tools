# action:
#     set environment variables for genome file access
#     index $GENOME.fa as needed
#     collect chromosome metadata
# expects:
#     $GENOME_DIR, will default to $MDI_DIR/resources/genomes/$GENOME
#     $GENOME
#     indexed fasta file $GENOME.fa or genome.fa in $GENOME_DIR
# usage:
#     source $MODULES_DIR/genome/set_genome_vars.sh

# file prefixes for genome-specific output files
export DATA_GENOME_PREFIX=${DATA_FILE_PREFIX}.${GENOME}

# set and check the genome directory
if [[ "$GENOME_DIR" == "" || "$GENOME_DIR" == "null" || "$GENOME_DIR" == "NA" ]]; then
    export GENOME_DIR=${MDI_DIR}/resources/genomes/${GENOME}
fi
if [ ! -d $GENOME_DIR ]; then
    echo
    echo "--genome-dir not found:"
    echo $GENOME_DIR
    echo
    exit 1
fi
export GENOME_PREFIX=${GENOME_DIR}/${GENOME}

# fasta file and index
export GENOME_FASTA=${GENOME_PREFIX}.fa
if [ ! -f $GENOME_FASTA ]; then
    export GENOME_FASTA=${GENOME_DIR}/genome.fa
fi
if [ ! -f $GENOME_FASTA ]; then
    echo "missing genome fasta file in ${GENOME_DIR}"
    echo "expected either file ${GENOME}.fa or genome.fa"
    exit 1
fi
if [ ! -f $GENOME_FASTA.fai ]; then
    echo "indexing genome fasta file"
    samtools index ${GENOME_FASTA}
fi

# metadata directories and files
export GENOME_METADATA_DIR=${GENOME_DIR}/metadata
export GENOME_METADATA_PREFIX=${GENOME_METADATA_DIR}/${GENOME}
export GENOME_GAPS_FILE=${GENOME_METADATA_PREFIX}.gaps.txt
export GENOME_EXCLUSIONS_BED=${GENOME_METADATA_PREFIX}.exclusions.bed
export GENOME_GC5BASE_WIG=${GENOME_METADATA_PREFIX}.gc5Base.wigVarStep.gz
export GENOME_BINS_DIR=${GENOME_DIR}/bins

# annotations
export GENOME_ANNOTATIONS_DIR=${GENOME_DIR}/annotations
export GENOME_GENCODE_PREFIX=${GENOME_ANNOTATIONS_DIR}/${GENOME}.gencode.${GENCODE_RELEASE}
export GENOME_GTF=${GENOME_GENCODE_PREFIX}.basic.annotation.gtf.gz
export GENES_BED=${GENOME_GENCODE_PREFIX}.genes.bed.gz

# get genome size from canonical chromosomes only
# get the list of all placed chromosome sequences, including chrX, chrY, chrM, and chrEBV if present
export GENOME_SIZE=`awk '$1!~/_/ && $1!="chrEBV"{s+=$2}END{print s}' $GENOME_FASTA.fai`
export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1 | grep -v _`
