#!/usr/bin/env python3 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr

# load in bedgraph files
def load_bedgraph(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])
    df['chrom'] = df['chrom'].astype(str)
    return df

# generate scatterplot of 2 bedgraph compartments
def plot_scatter(df1, df2, out_dir, res):

    plt.figure(figsize=(8, 8))
    plt.scatter(df1['value'], df2['value'], alpha=0.5)
    plt.xlabel('Bedgraph 1 Values')
    plt.ylabel('Bedgraph 2 Values')
    plt.title('Scatter Plot of Bedgraph Compartments')
    plt.grid()
    plt.show()
    # save plot to output dir

    plt.savefig(f'{out_dir}/scatter_plot_bedgraph_compartments_{str(res)}.png')

def main():
    out_dir = '/nfs/turbo/umms-hammou/minjilab/juicer/plots/121125_compartments_all_samples_combined_vs_ESE14'
    # Create the directory and any necessary parent directories
    # if they don't already exist.
    os.makedirs(out_dir, exist_ok=True)
    # Load bedgraph files
    #df2 = load_bedgraph("/nfs/turbo/umms-hammou/minjilab/juicer/work/mega_subset/aligned/compartments/flipped_compartments_250000_pc1.bedGraph")
    df1 = load_bedgraph("/nfs/turbo/umms-hammou/minjilab/juicer/work/ESE14/4DN_ESE14_compartments.bedGraph")
    #df1 = load_bedgraph("/nfs/turbo/umms-hammou/minjilab/juicer/work/ESE14/compartments/compartments_250000_pc1.bedGraph")
    #df2 = load_bedgraph("/nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/mega/aligned/compartments/inter_30_converted_hictk_AB_aligned.bedGraph")
    df2 = load_bedgraph("/nfs/turbo/umms-hammou/minjilab/juicer/work/all_samples_combined/mega/aligned/compartments/flipped_compartments_250000_pc1.bedGraph")
    df1 = df1[df1.chrom != 'chrY']
    df2 = df2[df2.chrom != 'chrY']
    #print(df1.head())
    #print(df2.head())
    #df1['start'] = df1['start'].astype(int)
    #df1['end'] = df1['end'].astype(int)

    #df2 = df2.reset_index()
    #df2.columns = ['chrom', 'start', 'end', 'value', 'comp']    
    #df2['start'] = df2['start'].astype(int)
    #df2['end'] = df2['end'].astype(int)

    # Merge on chrom, start, end
    merged = pd.merge(df1, df2, on=["chrom", "start", "end"], suffixes=("_1", "_2"))
    if merged.empty:
        print("No overlapping bins between the two bedGraph files.")
        return
    
    # get resolution from start-end of first row
    res = merged['end'][0] - merged['start'][0]
    
    merged['value_1'] = merged['value_1'].fillna(0)
    merged['value_2'] = merged['value_2'].fillna(0)
    
    # Normalize value_1 and value_2 to -1 to 1, keeping sign
    for col in ['value_1', 'value_2']:
        max_abs = merged[col].abs().max()
        if max_abs != 0:
            merged[col] = merged[col] / max_abs
        else:
            merged[col] = 0  # If all values are zero, set to 0
    

    # Compute correlation
    corr, _ = pearsonr(merged['value_1'], merged['value_2'])
    print(f"Correlation between bedgraph1 and bedgraph2: {corr}")

    # Plot values for all chromosomes present in merged dataframe
    chromosomes = merged['chrom'].unique()
    for chrom_to_plot in chromosomes:
        chrom_df = merged[merged['chrom'] == chrom_to_plot].copy()
        if chrom_df.empty:
            print(f"No data for {chrom_to_plot} in merged bedGraphs.")
            continue
        # Compute correlation for this chromosome
        chrom_corr, _ = pearsonr(chrom_df['value_1'], chrom_df['value_2'])
        print(f"Correlation for {chrom_to_plot}: {chrom_corr}")
        plt.figure(figsize=(16, 5))
        plt.plot(chrom_df['start'], chrom_df['value_1'], label='ESE14', alpha=0.7)
        plt.plot(chrom_df['start'], chrom_df['value_2'], label='Combined', alpha=0.7)
        plt.xlabel('Genomic position (start)')
        plt.ylabel('Compartment value')
        plt.title(f'{res} Resolution Compartment values for {chrom_to_plot} (corr={chrom_corr:.3f})')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{out_dir}/compartment_tracks_{chrom_to_plot}_{res}.png')
        plt.close()
    print(f"Saved compartment value track plots for chromosomes: {', '.join(chromosomes)}")


if __name__ == "__main__":
    main()