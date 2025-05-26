#!/bin/bash

# set insert size-dependent thresholds for number of errors
get_n_errors() {
    local kmer_length=$1
    if [ $kmer_length -le 40 ]; then
        echo 1
    elif [ $kmer_length -le 100 ]; then
        echo 2
    elif [ $kmer_length -le 200 ]; then
        echo 3
    else
        echo 4
    fi
}

# create the genome index once
if [ -f ${GENOME_GENMAP_DIR}/index/index.ids.concat ]; then
    echo "genmap index files already exist, skipping index creation"
else 
    echo "indexing $GENOME fasta file for genmap"
    rm -rf ${GENOME_GENMAP_DIR}/index
    genmap index -F ${GENOME_FASTA} -I ${GENOME_GENMAP_DIR}/index
    checkPipe
fi
echo

# create a mappability map for each kmer length
echo "creating mappability maps for $GENOME"
mkdir -p ${GENOME_GENMAP_DIR}/maps
for KMER_LENGTH in $MAPPABILITY_KMER_LENGTHS; do
    N_ERRORS=$(get_n_errors $KMER_LENGTH)
    echo "insert size = $KMER_LENGTH, allowed errors = $N_ERRORS"
    OUTPUT_PREFIX=${GENOME_GENMAP_DIR}/maps/${GENOME}.mappability.k_${KMER_LENGTH}.e_${N_ERRORS}
    export KMER_LENGTH=${KMER_LENGTH} # for use in perl script

    if [ -f ${OUTPUT_PREFIX}.bed.gz ]; then
        echo "    genmap bed file already exists, skipping kmer length"
    else 

        # run genmap
        if [ ! -f ${OUTPUT_PREFIX}.bedgraph ]; then
            genmap map \
                --threads $N_CPU \
                -K ${KMER_LENGTH} -E ${N_ERRORS} \
                -I ${GENOME_GENMAP_DIR}/index \
                -O ${OUTPUT_PREFIX} \
                -bg
            checkPipe
        fi

        # simplify to runs 1(true) = mappable
        # everything not in a run is implied 0(false) (unmappable)
        cat ${OUTPUT_PREFIX}.bedgraph |
        awk '$4==1' |

        # split runs at bin boundaries
        # extend the end of the run to include the bases to the end of the rightmost kmer
        perl $ACTION_DIR/split_at_bins.pl | 



        # TODO: consider here collecting sequences with bedtools getfasta (instead bedtools nuc)
        # and using it to calculate the GC content of the kmers AND collect the Tn5 kmers
        # thus, downstream code will be able to calculate Tn5 site weights per bin

        # or do this in bins.sh

        bedtools getfasta -bedOut -fi $FA_FILE -bed - | 
        perl $ACTION_DIR/parse_bin_sequences.pl | 



        # calculate the percent GC content of the bases contributing to each bin-split run
        bedtools nuc -fi ${GENOME_FASTA_SHM} -bed -  |
        tail -n +2 | # remove the bedtools nuc header

        # adjust the end positions back so that mappability runs only refer to left-aligned kmer starts
        # extract bedtools nuc GC column and round to 4 decimal places
        awk 'BEGIN{OFS="\t"}{ print $1, $2, $3 - '${KMER_LENGTH}' + 1, int($5 * 10000 + 0.5) / 10000 }' |

        # write the compressed bed file
        # contains runs of the left-aligned position of mappable kmers with fraction GC of the contributing bases
        # runs of mappable kmers are split at bin boundaries, but the bases contributing to a run may cross bin boundaries
        pigz -c --processes ${N_CPU} > ${OUTPUT_PREFIX}.bed.gz

        # # remove the uncompressed file
        # rm -f ${OUTPUT_PREFIX}.bedgraph
    fi
done

echo
echo "done"
