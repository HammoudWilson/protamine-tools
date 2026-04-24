#!/bin/bash
#SBATCH --job-name=merge_hic_juicer_all_samples
#SBATCH --account=hammou0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=12
#SBATCH --output=output/merge_hic_all_samples_%j.out
#SBATCH --error=error/merge_hic_all_samples_%j.err
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=160G
#SBATCH --time=10:00:00
#SBATCH --profile=Task

# >>> Conda initialization >>>
source ~/.bashrc
conda activate juicer_env
# <<< Conda initialization <<<


BASE_DIR=/nfs/turbo/umms-hammou/minjilab/juicer
GENOME=mm10
REFERENCE=$BASE_DIR/references/mm10_index/mm10_no_alt_analysis_set_ENCODE.fa
CHROMSIZES=$BASE_DIR/references/mm10_no_alt.chrom.sizes_ENCSR425FOI.tsv
SITES=$BASE_DIR/restriction_sites/mm10_GATC_GANTC.txt.gz


# Create a new folder for the mega run
mkdir -p /nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined

# Create symlinks to the desired directories
ln -s /nfs/turbo/umms-hammou/minjilab/juicer/work/day35 /nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/070325_day35
ln -s /nfs/turbo/umms-hammou/minjilab/juicer/work/day39 /nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/070325_day39
ln -s /nfs/turbo/umms-hammou/minjilab/juicer/work/day35_sample1 /nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/110625_day35_sample1
ln -s /nfs/turbo/umms-hammou/minjilab/juicer/work/day35_sample2 /nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/110625_day35_sample2



# Run mega.sh to merge all merged_nodups.txt files and create .hic file
TOPDIR="$BASE_DIR/work/all_samples_combined"  # Set to your top-level directory


cd $TOPDIR
bash /nfs/turbo/umms-hammou/minjilab/juicer/scripts/mega.sh -g $GENOME -s $SITES \
    -D $BASE_DIR \
    -q standard \
    -l largemem \
    -L "20:00:00" \
    -Q "10:00:00" 

echo "mega.sh run completed. Merged files and .hic file should be in $TOPDIR/mega/aligned/"
