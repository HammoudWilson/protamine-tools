#!/bin/bash

# set the composite chrom name delimiter
DELIMITER="-"

# composite genome FASTA
echo
echo "concatenating genome FASTA files"
cat <(
    cat ${PRIMARY_GENOME_FASTA} | 
    awk '{if($1 ~ /^>/) { print $1"'${DELIMITER}${PRIMARY_GENOME}'" }  else { print $0} }'
) <(
    cat ${SPIKE_IN_GENOME_FASTA} | 
    awk '{if($1 ~ /^>/) { print $1"'${DELIMITER}${SPIKE_IN_GENOME}'" } else { print $0} }'
) > ${GENOME_FASTA}
checkPipe

# composite genome FASTA index
echo "indexing concatenated genome FASTA file"
samtools faidx ${GENOME_FASTA}
checkPipe

# composite gaps file
echo "concatenating genome gaps files"
cat <(
    awk 'BEGIN{OFS="\t"}{$2 = $2"'${DELIMITER}${PRIMARY_GENOME}'";  print $0}' ${PRIMARY_GENOME_GAPS_FILE}
) <(
    awk 'BEGIN{OFS="\t"}{$2 = $2"'${DELIMITER}${SPIKE_IN_GENOME}'"; print $0}' ${SPIKE_IN_GENOME_GAPS_FILE}
) > ${GENOME_GAPS_FILE}
checkPipe

# composite exclusions file
echo "concatenating genome exclusion files"
cat <(
    cat ${PRIMARY_GENOME_EXCLUSIONS_BED} |
    cut -f 1-3 |
    awk 'BEGIN{OFS="\t"}{$1 = $1"'${DELIMITER}${PRIMARY_GENOME}'";  print $0}'
) <(
    cat ${SPIKE_IN_GENOME_EXCLUSIONS_BED} |
    cut -f 1-3 |
    awk 'BEGIN{OFS="\t"}{$1 = $1"'${DELIMITER}${SPIKE_IN_GENOME}'"; print $0}'
) > ${GENOME_EXCLUSIONS_BED}
checkPipe

# composite inclusions file
echo "parsing gaps and exclusions files into an inclusion regions file"
cat <(cat ${GENOME_GAPS_FILE} | cut -f 2-4) ${GENOME_EXCLUSIONS_BED} |
sort -k1,1 -k2,2n -k3,3n |
bedtools complement -i stdin -g <(cat ${GENOME_FASTA}.fai | sort -k1,1) |
awk '$3 - $2 >= 2000' | # guarantees at least one 1kb bin per inclusion region
grep -v "_" > ${GENOME_INCLUSIONS_BED}
checkPipe

# composite genes bed (primary genome genes only)
echo "tranferring primary genome annotation BED"
zcat ${PRIMARY_GENOME_GENES_BED} |
awk 'BEGIN{OFS="\t"}{$1 = $1"'${DELIMITER}${PRIMARY_GENOME}'"; print $0}' |
gzip -c > ${GENES_BED}
checkPipe

echo
echo "done"
echo
