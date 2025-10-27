
# set any constants/environment variables needed by the script
ROUND_M=10000 # multiplier when rounding percent GC values

# declare the input file(s) used to create the output BED file
# this is the baseline genome 1kb bins file created by previous genome/bin with single bin GC value
#    chrom   start0  end1    excluded        pct_gc
EXT_GENOMES_DIR=/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/resources/genomes
INPUT_FILE=${EXT_GENOMES_DIR}/${PRIMARY_GENOME}/bins/${PRIMARY_GENOME}.bin_1000.bed.gz

# name and excerpt the input file(s) to the job log
# be sure to include the redirection to STDERR (1>&2) on all log file lines!
# choose `head` or `zcat` depending on whether the input file is compressed
echo "input file: ${INPUT_FILE}" 1>&2
zcat ${INPUT_FILE} | head 1>&2

# keep all bins as is with no merging, round pct_gc values
# this provides the most accurate local value for a small query region at the expense of a large file
zcat ${INPUT_FILE} |
awk '$1 != "chrom"' | # remove the header
awk -v M=${ROUND_M} 'BEGIN{
    OFS = "\t";
    print "chrom", "start0", "end1", "name", "Fraction_GC";
} {
    frac_gc = int($5 * M + 0.5) / M;
    print $1, $2, $3, ".", frac_gc;
}'
