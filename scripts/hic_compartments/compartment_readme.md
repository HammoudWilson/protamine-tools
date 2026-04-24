Merge all of the hic directories together using merge_hic_juicer.sh

Command line order of running to flip the compartments:
sbatch compartment_calling.sh -o /nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/mega/aligned/compartments -h work/all_samples_combined/mega/aligned/inter_30.hic
./flip_compartments.py
./compartment_correlation.py

Liftover from mm10 to mm39 using liftover_mm10_to_mm39.sh