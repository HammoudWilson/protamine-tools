# declare the input file(s) used to create the output BED file
# this is the hic compartments called on merged early and late round spermatids at 1mb binsize
#    chrom   start0  end1    excluded        compartment_score
INPUT_FILE="/nfs/turbo/umms-hammou/minjilab/juicer/work/mega_subset/aligned/compartments/compartments_1000000_pc1.bedGraph"

# name and excerpt the input file(s) to the job log
# be sure to include the redirection to STDERR (1>&2) on all log file lines!
# choose `head` or `zcat` depending on whether the input file is compressed
echo "input file: ${INPUT_FILE}" 1>&2
head -n 2 "${INPUT_FILE}" 1>&2

# Add header and output as BED format 
sed "1i chrom\tstart0\tend1\tcompartment_score" "${INPUT_FILE}"
