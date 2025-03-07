#!/bin/bash

# This script concatenates multi-lane FastQ files for each sample into a single R1 and R2 file. 
# It takes raw FastQ files from a specified directory, merges lanes, and outputs them to a processed directory.
# Sample names are mapped to standardized identifiers to ensure consistent naming in downstream analyses.

# Define input and output directories
SAMPLEDIR="/path/to/raw_fastq"  # Directory containing raw FastQ files
ENDDIR="/path/to/processed_fastq"  # Directory where merged FastQ files will be stored

cd "$SAMPLEDIR"

# Define original sample names and their corresponding standardized identifiers
declare -a sample_names=("Sample1_S1" "Sample2_S2" "Sample3_S3" "Sample4_S4")
declare -a sample_rename=("sampleA" "sampleB" "sampleC" "sampleD")

# Ensure the output directory exists
mkdir -p "$ENDDIR"

# Loop through each sample, concatenate lane-specific FastQ files, and save them with standardized names
echo "Processing samples..."
for (( i=0; i<${#sample_names[@]}; i++ )); do
    echo "Merging lanes for: ${sample_rename[i]}"
    
    cat "${sample_names[i]}_L001_R1_001.fastq.gz" \
        "${sample_names[i]}_L002_R1_001.fastq.gz" \
        "${sample_names[i]}_L003_R1_001.fastq.gz" \
        "${sample_names[i]}_L004_R1_001.fastq.gz" > "$ENDDIR/${sample_rename[i]}_R1.fastq.gz"

    cat "${sample_names[i]}_L001_R2_001.fastq.gz" \
        "${sample_names[i]}_L002_R2_001.fastq.gz" \
        "${sample_names[i]}_L003_R2_001.fastq.gz" \
        "${sample_names[i]}_L004_R2_001.fastq.gz" > "$ENDDIR/${sample_rename[i]}_R2.fastq.gz"

    echo "Merged ${sample_rename[i]} R1 and R2 FastQ files."
done

echo "FastQ merging complete."

