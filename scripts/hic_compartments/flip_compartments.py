#!/usr/bin/env python3 
import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr

align_track = pd.read_csv("/nfs/turbo/umms-hammou/minjilab/juicer/work/ESE14/4DNFIPUC6WYV_ESE14_compartments.bedGraph", sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])

#compartment_files = ["/nfs/turbo/umms-hammou/minjilab/juicer/work/mega_subset/aligned/compartments/compartments_250000_pc1.bedGraph",
#                     "/nfs/turbo/umms-hammou/minjilab/juicer/work/day35/aligned/compartments/compartments_250000_pc1.bedGraph",
#                     "/nfs/turbo/umms-hammou/minjilab/juicer/work/day39/aligned/compartments/compartments_250000_pc1.bedGraph"]

compartment_files = ["/nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/mega/aligned/compartments/compartments_250000_pc1.bedGraph"]

# add chr to first column of align_track
#align_track['chrom'] = 'chr' + align_track['chrom'].astype(str)

for compartment_file in compartment_files:
    compartment_track = pd.read_csv(compartment_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])


    # group by chromosome and calculate correlation
    grouped = compartment_track.groupby('chrom')
    df1 = align_track[align_track.chrom != 'chrY']
    df2 = compartment_track[compartment_track.chrom != 'chrY']

    # Merge on chrom, start, end
    merged = pd.merge(df1, df2, on=["chrom", "start", "end"], suffixes=("_1", "_2"))
    merged = merged.dropna()
    if merged.empty:
        print("No overlapping bins between the two bedGraph files.")
        continue
    #print(merged.head())
    # Compute correlation for each chromosome
    chroms = merged['chrom'].unique()
    for chrom in chroms:
        chrom_data = merged[merged['chrom'] == chrom]
        #print(chrom_data.head())
        corr, _ = pearsonr(chrom_data['value_1'], chrom_data['value_2'])
        #print(f"Chromosome: {chrom}, Correlation: {corr}")

        if corr < 0:
            # Flip the sign of the compartment values for that chromosome
            merged.loc[merged['chrom'] == chrom, 'value_2'] *= -1

    # get chrom, start, end, and value_2 columns for bedgraph output
    output = merged[['chrom', 'start', 'end', 'value_2']]
    output_file = os.path.join(os.path.dirname(compartment_file), f"flipped_{os.path.basename(compartment_file)}")
    output.to_csv(output_file, sep='\t', header=False, index=False)

    
    
    
    