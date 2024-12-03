#!/bin/bash

# Loop through all BAM files ending with .mm39.dedup.bam in the input directory
for bam_file in "$INPUT_DIR"/*.mm39.dedup.bam; do

    # Extract the base name of the BAM file (without the directory and extension)
    base_name=$(basename "$bam_file" .mm39.dedup.bam)

    # Define the output file path
    output_file="$OUTPUT_DIR/${base_name}_bam_direct_cov.bed"

    # Run bedtools genomecov and save the output
    bedtools coverage -counts -a ${GENOME_BINS_BED} -b "$bam_file" -sorted > "$output_file"

    echo "Processed $bam_file -> $output_file"
done


# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("~/turbo_hammoud/DATA_STORAGE/ATAC-processed/Tim_Parnell_analysis/A8378_10%_drosophila_spikein/mashiat_processed/bam_direct_intersect_windows/")

# Load filenames matching the pattern
filenames <- list.files(path = '.', pattern = 'bam_direct_cov_500kb.bed', full.names = TRUE)
filenames <- sub("^\\./", "", filenames) # Remove './' from filenames
print(filenames)

# Function to read and format a single file
read_data <- function(filepath) {
  file <- read.delim2(filepath, header = FALSE)
  coverage_name <- sub("_bam_direct_cov_500kb.bed$", "", basename(filepath))
  coverage_name <- paste0('s_', coverage_name)
  colnames(file) <- c('chr', 'start', 'end', coverage_name)
  return(file)
}

# Read all files into a list of data frames
coverage_list <- lapply(filenames, read_data)

# Dynamically rename each data frame using a sanitized version of the filename
names(coverage_list) <- sapply(filenames, function(filepath) {
  sub("_bam_direct_cov_500kb.bed$", "", basename(filepath))
})

# Combine all data frames into a single data frame by full joining
coverage_combined <- Reduce(function(x, y) full_join(x, y, by = c("chr", "start", "end")), coverage_list)

# Print the head of the combined data frame
print(head(coverage_combined))

# Write the combined data frame to a CSV file
write.csv(coverage_combined, file = 'bam_coverage_500kb_combined.csv', row.names = FALSE)





