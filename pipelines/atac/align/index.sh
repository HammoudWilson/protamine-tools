# action:
#     create the Bowtie2 index of the composite genome if needed
# input:
#     composite genome fasta file
# outputs:
#     bowtie2 index files

if [ -f ${GENOME_BOWTIE2_PREFIX}.1.bt2 ]; then
    echo "bowtie2 index files already exist, skipping index creation"
else 
    echo "creating bowtie2 index files for $GENOME"
    bowtie2-build --threads ${N_CPU} --quiet ${GENOME_FASTA} ${GENOME_BOWTIE2_PREFIX}
    checkPipe
    echo "done"
fi
