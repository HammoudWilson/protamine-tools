#!/bin/bash

echo "creating genome bins..."

# create fixed-width genome bins
bedtools makewindows -g ${GENOME_FASTA}.fai -w ${BIN_SIZE} |

# mark (but do not delete) bins that are excluded from analysis as (partially) overlapping with:
#   genome gaps
#   problematic regions
bedtools intersect -c -a - -b ${GENOME_EXCLUSIONS_BED} <(cut -f2-4 ${GENOME_GAPS_FILE}) | 
awk '{print $0"\t"($NF>0 ? 1 : 0)}' | 
cut -f1-3,5 |

# calculate the percent GC of each bin
bedtools nuc -fi ${GENOME_FASTA} -bed - |
awk 'NR>1' |
cut -f1-4,6 |

# sort bins in a manner consistent with downstream steps
sort --parallel ${N_CPU} --buffer-size 1G -k1,1V -k2,2n |

# add a header and save the final gzipped file
awk 'BEGIN{print "chrom\tstart\tend\texcluded\tpct_gc"} 1' |
pigz -c --processes ${N_CPU} > ${GENOME_BINS_BED}
checkPipe

echo "done"
echo
echo "head of ${GENOME_BINS_BED}:"
zcat ${GENOME_BINS_BED} | head
