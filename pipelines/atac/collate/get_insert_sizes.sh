# action:
#     calculate insert size distributions in a single bam file
# input:
#     $BAM_FILE as arg1
#     $GENOME_FASTA_SHM, on /dev/shm for maxmimum lookup performance
#     $GENOME_INCLUSIONS_BED
#     $MIN_MAPQ
#     $MIN_INSERT_SIZE
#     $MAX_INSERT_SIZE
#     $PRIMARY_GENOME
#     $SPIKE_IN_GENOME
# outputs:
#     sample-level insert size distribution table on STDOUT
#       rows    = insert sizes
#       columns = genome, percent GC rounded from 0 to 100

# get arguments
export BAM_FILE=$1
export CHROM="" # thus, this script acts over all chromosomes

# pull all insert endpoints from the bam file
# since final output is aggregated, enforce inclusions here rather than downstream
bash $ACTION_DIR/../pull_insert_endpoints.sh ENFORCE_INCLUSIONS | 

# collect the GC content of the inserts
# this is done from the reference genome, not the read, to ensure 
#   that GC content of larger paired-read inserts reflects their entire span, not just the reads
bedtools nuc -fi $GENOME_FASTA_SHM -bed - |
tail -n +2 | # remove the bedtools nuc header
cut -f 1-3,5 |

# ###########################
# # head -n 5000000 | 

# deduplicate and count insert sizes stratified by GC content
perl $ACTION_DIR/get_insert_sizes.pl
