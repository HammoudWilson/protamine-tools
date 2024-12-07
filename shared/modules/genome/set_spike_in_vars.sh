# action:
#     set environment variables for spike-in file access
#     index $SPIKE_IN_GENOME as needed
#     collect chromosome metadata
# expects:
#     $SPIKE_IN_GENOME_DIR, will default to $MDI_DIR/resources/genomes/$SPIKE_IN_GENOME
#     $SPIKE_IN_GENOME
#     indexed fasta file $SPIKE_IN_GENOME.fa or genome.fa in $SPIKE_IN_GENOME_DIR
# usage:
#     source $MODULES_DIR/genome/set_spike_in_vars.sh

# file prefixes for spike-in-specific output files
export DATA_SPIKE_IN_PREFIX=${DATA_FILE_PREFIX}.${SPIKE_IN_GENOME}

# set and check the spike-in directory
if [[ "$SPIKE_IN_GENOME_DIR" == "" || "$SPIKE_IN_GENOME_DIR" == "null" || "$SPIKE_IN_GENOME_DIR" == "NA" ]]; then
    export SPIKE_IN_GENOME_DIR=${MDI_DIR}/resources/genomes/${SPIKE_IN_GENOME}
fi
if [ ! -d $SPIKE_IN_GENOME_DIR ]; then
    echo
    echo "--spike-in-genome-dir not found:"
    echo $SPIKE_IN_GENOME_DIR
    echo
    exit 1
fi
export SPIKE_IN_PREFIX=${SPIKE_IN_GENOME_DIR}/${SPIKE_IN_GENOME}

# fasta file and index
export SPIKE_IN_FASTA=${SPIKE_IN_PREFIX}.fa
if [ ! -f $SPIKE_IN_FASTA ]; then
    export SPIKE_IN_FASTA=${SPIKE_IN_GENOME_DIR}/genome.fa
fi
if [ ! -f $SPIKE_IN_FASTA ]; then
    echo "missing spike-in-genome fasta file in ${SPIKE_IN_GENOME_DIR}"
    echo "expected either file ${SPIKE_IN_GENOME}.fa or genome.fa"
    exit 1
fi
if [ ! -f $SPIKE_IN_FASTA.fai ]; then
    echo "indexing spike-in-genome fasta file"
    samtools index ${SPIKE_IN_FASTA}
fi

# metadata directories and files
export SPIKE_IN_METADATA_DIR=${SPIKE_IN_GENOME_DIR}/metadata
export SPIKE_IN_METADATA_PREFIX=${SPIKE_IN_METADATA_DIR}/${SPIKE_IN_GENOME}
export SPIKE_IN_GAPS_FILE=${SPIKE_IN_METADATA_PREFIX}.gaps.txt
export SPIKE_IN_EXCLUSIONS_BED=${SPIKE_IN_METADATA_PREFIX}.exclusions.bed
export SPIKE_IN_GC5BASE_WIG=${SPIKE_IN_METADATA_PREFIX}.gc5Base.wigVarStep.gz
export SPIKE_IN_BINS_DIR=${SPIKE_IN_GENOME_DIR}/bins

# annotations
export SPIKE_IN_ANNOTATIONS_DIR=${SPIKE_IN_GENOME_DIR}/annotations

# get spike-in-genome size from canonical chromosomes only
# get the list of all placed chromosome sequences, including chrX, chrY, chrM, and chrEBV if present
export SPIKE_IN_SIZE=`awk '$1!~/_/ && $1!="chrEBV"{s+=$2}END{print s}' $SPIKE_IN_FASTA.fai`
export SPIKE_IN_CHROMS=`cat $SPIKE_IN_FASTA.fai | cut -f1 | grep -v _`
