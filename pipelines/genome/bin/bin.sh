#!/bin/bash

echo "creating genome bins based on:"
echo "  ${GENOME_FASTA}"
echo "  ${GENOME_GAPS_FILE}"
echo "  ${GENOME_EXCLUSIONS_BED}"
echo "  bin size = ${BIN_SIZE}"

# filter and sort chromosomes in a manner consistent with downstream steps
cut -f 1-2 ${GENOME_FASTA}.fai | 
grep -v -e "_" -e "chrM" -e "chrEBV" | 
sort -k1,1V |

# create fixed-width genome bins
bedtools makewindows -g - -w ${BIN_SIZE} |

# mark (but do not delete) bins that are excluded from analysis as (partially) overlapping with:
#   genome gaps
#   problematic regions
bedtools intersect -c -a - -b ${GENOME_EXCLUSIONS_BED} <(cut -f2-4 ${GENOME_GAPS_FILE}) | 
awk '{print $0"\t"($NF>0 ? 1 : 0)}' | # thus, output is: chrom,start0,end1,excluded=[0,1]
cut -f1-3,5 |

# calculate the percent GC of each bin
# bedtools nuc output is original bed columns (four from above) plus:
#     %AT content (pct_at)
#     %GC content (pct_gc) (actually expressed as a fraction)
bedtools nuc -fi ${GENOME_FASTA} -bed - |
awk 'NR>1' |
cut -f1-4,6 |

# add a header and save the final gzipped file
awk 'BEGIN{print "chrom\tstart0\tend1\texcluded\tpct_gc"} 1' |
pigz -c --processes ${N_CPU} > ${GENOME_BINS_BED}
checkPipe

echo "done"
echo
echo "head of ${GENOME_BINS_BED}:"
zcat ${GENOME_BINS_BED} | head
echo
echo "chrom ranges"
zcat ${GENOME_BINS_BED} | 
bedtools groupby -g 1 -c 2,3,4 -o first,last,mean
