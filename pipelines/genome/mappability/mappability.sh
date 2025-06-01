# action:
#     use genmap to assess genome mappability at different insert size levels
#     mappability can be different (lower) for smaller insert sizes seen in ATAC-seq data
# input:
#     $GENOME_FASTA
# output:
#     a series of compressed bed files with genome mappability run spans and GC content
#       one file per insert size level
#       chrom, start0_run, end1_run, fraction_GC
#       where _run indicates the start and end of a run of mappable insert starts

# set isert size-dependent thresholds for number of errors
get_n_errors() {
    local INSERT_SIZE_LEVEL=$1
    if [ ${INSERT_SIZE_LEVEL} -le 40 ]; then
        echo 1
    elif [ ${INSERT_SIZE_LEVEL} -le 100 ]; then
        echo 2
    elif [ ${INSERT_SIZE_LEVEL} -le 200 ]; then
        echo 3
    else
        echo 4
    fi
}

# create the genome index once
if [ -f ${GENOME_GENMAP_DIR}/index/index.ids.concat ]; then
    echo "genmap index files already exist, skipping index creation"
else 
    echo "indexing ${GENOME} fasta file for genmap"
    rm -rf ${GENOME_GENMAP_DIR}/index
    genmap index -F ${GENOME_FASTA} -I ${GENOME_GENMAP_DIR}/index
    checkPipe
fi
echo

# create a mappability map for each insert size level
echo "creating mappability maps for ${GENOME}"
mkdir -p ${GENOME_GENMAP_DIR}/maps
for INSERT_SIZE_LEVEL in ${MAPPABILITY_SIZE_LEVELS}; do
    N_ERRORS=$(get_n_errors ${INSERT_SIZE_LEVEL})
    echo "insert size level = ${INSERT_SIZE_LEVEL}, allowed errors = ${N_ERRORS}"
    OUTPUT_PREFIX=${GENOME_GENMAP_DIR}/maps/${GENOME}.mappability.k_${INSERT_SIZE_LEVEL}.e_${N_ERRORS}
    export INSERT_SIZE_LEVEL=${INSERT_SIZE_LEVEL} # for use in perl script

    if [ -f ${OUTPUT_PREFIX}.bed.gz ]; then
        echo "    genmap bed file already exists, skipping insert size level"
    else 

        # run genmap
        if [ ! -f ${OUTPUT_PREFIX}.bedgraph ]; then
            genmap map \
                --threads ${N_CPU} \
                -K ${INSERT_SIZE_LEVEL} -E ${N_ERRORS} \
                -I ${GENOME_GENMAP_DIR}/index \
                -O ${OUTPUT_PREFIX} \
                -bg
            checkPipe
        fi

        # simplify to runs 1(true) = mappable
        # everything not in a run is implied 0(false) = unmappable
        cat ${OUTPUT_PREFIX}.bedgraph |
        awk '$4==1' |

        # split runs at bin boundaries
        # extend the end of each run to include the bases to the end of the rightmost possible insert
        perl ${ACTION_DIR}/split_at_bins.pl | 

        # calculate the percent GC content of the bases contributing to each bin-split run
        bedtools nuc -fi ${GENOME_FASTA_SHM} -bed -  |
        tail -n +2 | # remove the bedtools nuc header

        # adjust the coordinate positions back so that mappability runs only refer to left-aligned insert starts
        # extract bedtools nuc GC column and round to 4 decimal places
        awk -v isl=${INSERT_SIZE_LEVEL} 'BEGIN{OFS="\t"}{ 
            print $1, $2, $3 - isl + 1, int($5 * 10000 + 0.5) / 10000 
        }' |

        # write the compressed bed file
        # contains runs of the left-aligned position of mappable inserts with fraction GC of the contributing bases
        # runs of mappable start positions are split at bin boundaries, but the bases contributing to a run may cross bin boundaries
        pigz -c --processes ${N_CPU} > ${OUTPUT_PREFIX}.bed.gz
        checkPipe

        # # remove the uncompressed file
        # rm -f ${OUTPUT_PREFIX}.bedgraph
    fi
done

echo
echo "done"
