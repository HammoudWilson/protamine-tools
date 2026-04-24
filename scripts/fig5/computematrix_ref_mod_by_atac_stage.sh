#!/bin/bash

#SBATCH --job-name=compute_atac_dinuc_early

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=4

#SBATCH --mem-per-cpu=32GB 

#SBATCH --time=40:00:00

#SBATCH --account=junzli1

#SBATCH --partition=standard

module load Bioinformatics

module load python/3.11.5



computeMatrix reference-point -S /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/Cut_Tag/A8672/26615R_histone/Fragment/day56_RS_H4.fragment.bw /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/Cut_Tag/A8672/26615R_histone/Fragment/day56_ES_H4.fragment.bw /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/Cut_Tag/A8672/26615R_histone/Fragment/day56_RS_H2b.fragment.bw /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/Cut_Tag/A8672/26615R_histone/Fragment/day56_ES_H2b.fragment.bw /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/Cut_Tag/A8538/25950R_histone/Fragment/D34_RS_H4ac.fragment.bw /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/Cut_Tag/A8538/25950R_histone/Fragment/D34_ES_H4ac.fragment.bw  --samplesLabel RS_H4 ES_H4 RS_H2B ES_H2B RS_H4ac ES_H4ac -R /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/ATAC-processed/Tim_Parnell_analysis/Mashiat_analysis/ATAC_alldatasets/bed_by_stage/dinuc_stage_Early.bed /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/ATAC-processed/Tim_Parnell_analysis/Mashiat_analysis/ATAC_alldatasets/bed_by_stage/dinuc_stage_mid.bed /nfs/turbo/umms-hammou/mrabbani/DATA_STORAGE/ATAC-processed/Tim_Parnell_analysis/Mashiat_analysis/ATAC_alldatasets/bed_by_stage/dinuc_stage_late.bed --referencePoint center -b 2500 -a 2500 --skipZeros --missingDataAsZero -o dinuc_atac_E_M_L_canonical_histones_h4_ac_rs_es.gz -p 8 --outFileNameMatrix dinuc_atac_E_M_L_canonical_histones_h4_ac_rs_es.tab --outFileSortedRegions dinuc_atac_E_M_L_canonical_histones_h4_ac_rs_es.bed

