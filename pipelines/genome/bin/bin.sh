#!/bin/bash

echo "creating genome bins based on:"
echo "  ${GENOME_FASTA}"
echo "  ${GENOME_GAPS_FILE}"
echo "  ${GENOME_EXCLUSIONS_BED}"
echo "  ${MAPPABILITY_FILE_PREFIX}.*"
echo "  bin size = ${BIN_SIZE}"

# create a temporary file for the bin data to which map_N and gc_N will be added
echo "calculating bins, inclusions, and GC content"
BIN_TMP_FILE=$(mktemp)

# filter chromosomes in a manner consistent with downstream steps
cut -f 1-2 ${GENOME_FASTA}.fai | 
grep -v -e "_" -e "chrM" -e "chrEBV" | 

# create fixed-width genome bins
bedtools makewindows -g - -w ${BIN_SIZE} |

# mark (but do not delete) bins that are included for analysis as not even partially overlapping with:
#   genome gaps
#   problematic regions
bedtools intersect -c -a - -b ${GENOME_EXCLUSIONS_BED} <(cut -f2-4 ${GENOME_GAPS_FILE}) | 
awk '{print $0"\t"($NF>0 ? 0 : 1)}' | # thus, output is: chrom,start0,end1,included=[0,1]
cut -f1-3,5 > ${BIN_TMP_FILE}

# for each k-mer length, for each bin, calculate
#   the fraction of k-mers left-aligned to bin that are mappable 
#   the aggregate GC content of those mappable k-mers, weighted over mappable run lengths
echo "appending mappability levels and GC content"
for KMER_LENGTH in ${MAPPABILITY_KMER_LENGTHS}; do
    echo "  insert size level = ${KMER_LENGTH}"
    MAP_N_FILE_GZ=$( ls -1 ${MAPPABILITY_FILE_PREFIX}.k_${KMER_LENGTH}.e_*.bed.gz )
    MAP_TMP_FILE=$(mktemp)
    paste ${BIN_TMP_FILE} <(
        bedtools intersect -wao -sorted \
            -a <(cut -f 1-3 ${BIN_TMP_FILE}) \
            -b <(zcat ${MAP_N_FILE_GZ} | grep -v -e "_" -e "chrM" -e "chrEBV") |
        perl ${ACTION_DIR}/calculate_bin_metrics.pl
    ) > ${MAP_TMP_FILE}
    mv ${MAP_TMP_FILE} ${BIN_TMP_FILE}
done

# add header with appropriate column names
HEADER="chrom\tstart0\tend1\tincluded"
for KMER_LENGTH in ${MAPPABILITY_KMER_LENGTHS}; do
    HEADER="${HEADER}\tmap_${KMER_LENGTH}\tgc_${KMER_LENGTH}"
done

# create final output
(echo -e "${HEADER}"; cat ${BIN_TMP_FILE}) |
pigz -c --processes ${N_CPU} > ${GENOME_BINS_BED}
checkPipe

# clean up
rm -f ${BIN_TMP_FILE}

# modify the summary sections to handle multiple k-mer lengths
echo
echo ${GENOME_BINS_BED}
echo
echo "head of bins file"
zcat ${GENOME_BINS_BED} | 
head

echo
echo "head of bins file (included regions only)"
echo -e "${HEADER}"
zcat ${GENOME_BINS_BED} | 
awk '$4==1' |
head

# create mappability header with appropriate column names
HEADER="mappability"
for KMER_LENGTH in ${MAPPABILITY_KMER_LENGTHS}; do
    HEADER="${HEADER}\tk_${KMER_LENGTH}"
done

# create table of mappability values (0 to 1, step 0.05)
echo
echo "mappability summary (included bins only)"
MAP_TMP_FILE=$(mktemp)
seq 0 0.05 1 > ${MAP_TMP_FILE}
for i in $(seq 1 20); do
    LEVEL_TMP_FILE=$(mktemp)
    paste ${MAP_TMP_FILE} <(
        zcat ${GENOME_BINS_BED} | 
        grep -v "chrom" |
        awk '$4==1' |
        awk '{print int($'$((5 + 2*(i-1)))' * 20 + 0.5) / 20}' |
        sort -k1,1n | 
        bedtools groupby -g 1 -c 1 -o count | 
        cut -f 2
    ) > ${LEVEL_TMP_FILE}
    mv ${LEVEL_TMP_FILE} ${MAP_TMP_FILE}
done
echo -e "${HEADER}"
cat ${MAP_TMP_FILE} 
rm -f ${MAP_TMP_FILE}

# create GC header with appropriate column names
HEADER="frac_gc"
for KMER_LENGTH in ${MAPPABILITY_KMER_LENGTHS}; do
    HEADER="${HEADER}\tk_${KMER_LENGTH}"
done

# create table of frac_gc values (0 to 1, step 0.05)
echo
echo "gc summary (included bins only)"
GC_TMP_FILE=$(mktemp)
seq 0 0.05 1 > ${GC_TMP_FILE}
for i in $(seq 1 20); do
    LEVEL_TMP_FILE=$(mktemp)
    paste ${GC_TMP_FILE} <(
        zcat ${GENOME_BINS_BED} | 
        grep -v "chrom" |
        awk '$4==1' |
        awk '{print int($'$((6 + 2*(i-1)))' * 20 + 0.5) / 20}' |
        sort -k1,1n | 
        bedtools groupby -g 1 -c 1 -o count | 
        cut -f 2
    ) > ${LEVEL_TMP_FILE}
    mv ${LEVEL_TMP_FILE} ${GC_TMP_FILE}
done
echo -e "${HEADER}"
cat ${GC_TMP_FILE} 
rm -f ${GC_TMP_FILE}


# # prepare properly sorted and filtered versions of the map_N files
# echo "presorting the mappability files"
# for KMER_LENGTH in ${MAPPABILITY_KMER_LENGTHS}; do
#     echo "  ${KMER_LENGTH}"
#     MAP_N_FILE_BGZ=$( ls -1 ${MAPPABILITY_FILE_PREFIX}.k_${KMER_LENGTH}.e_*.bed.bgz )
#     MAP_N_FILE_GZ=${MAPPABILITY_FILE_PREFIX}.k_${KMER_LENGTH}.bin_sorted.gz
#     if [ ! -f ${MAP_N_FILE_GZ} ]; then
#         zcat ${MAP_N_FILE_BGZ} | 
#         grep -v -e "_" -e "chrM" -e "chrEBV" | 
#         sort --parallel ${N_CPU} --buffer-size 4G -k1,1 -k2,2n | # same sort as below
#         pigz -c --processes ${N_CPU} > ${MAP_N_FILE_GZ}
#         checkPipe
#     fi
# done

# # filter and sort chromosomes in a manner consistent with downstream steps
# # this lexiocgraphic sort intermingles primary and spike-in chromosomes
# cut -f 1-2 ${GENOME_FASTA}.fai | 
# grep -v -e "_" -e "chrM" -e "chrEBV" | 
# sort -k1,1 |

# # calculate the percent GC of each bin
# # bedtools nuc output is original bed columns (four from above) plus:
# #     %AT content (pct_at)
# #     %GC content (pct_gc) (actually expressed as a fraction)
# bedtools nuc -fi ${GENOME_FASTA} -bed - |
# awk 'NR>1' | # remove the bedtools nuc header
# cut -f1-4,6 |
# awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,sprintf("%.3f", $5)}' > ${BIN_TMP_FILE}  # round percent GC to 3 decimal places

# # bin data through GC content now held in BIN_TMP_FILE

# # for each k-mer length, for each bin, calculate
# #   the fraction of k-mers that are mappable 
# #   the weighted average of the GC content of mappable k-mers
# echo "appending mappability levels and GC content"
# for KMER_LENGTH in ${MAPPABILITY_KMER_LENGTHS}; do
#     MAP_N_FILE_GZ=${MAPPABILITY_FILE_PREFIX}.k_${KMER_LENGTH}.bin_sorted.gz
#     MAP_TMP_FILE=$(mktemp)
#     paste ${BIN_TMP_FILE} <(
#         bedtools intersect -wao -sorted -a <(cut -f 1-3 ${BIN_TMP_FILE}) -b ${MAP_N_FILE_GZ} |
#         bedtools groupby -g 1,2,3 -c 7 -o sum |
#         awk '{print sprintf("%.3f", $4/($3-$2))}'
#     ) > ${MAP_TMP_FILE}
#     mv ${MAP_TMP_FILE} ${BIN_TMP_FILE}
# done
