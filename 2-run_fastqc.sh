#!/bin/bash
#SBATCH --job-name=fastqc_analysis
#SBATCH --mem=8G
#SBATCH --array=13-13%1

# This script runs FastQC quality control analysis on paired-end FASTQ files.
# It retrieves sample names from a text file and processes each sample sequentially using SLURM job arrays.
# The results are saved in a specified output directory.

# Load FastQC module if required (Uncomment if using a module system)
# module load fastqc

# Define paths
INPUT_DIR="/path/to/processed_fastq"
WORKDIR="/path/to/analysis_directory"
SAMPLE_FILE="$WORKDIR/sample_names.txt"
OUTPUT_DIR="$WORKDIR/pre_qc_reports"

# Read the sample name from the file based on the SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_FILE")

# Define input FASTQ files
FASTQ_R1="${INPUT_DIR}/${SAMPLE}_R1.fastq.gz"
FASTQ_R2="${INPUT_DIR}/${SAMPLE}_R2.fastq.gz"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Run FastQC
echo "Running FastQC on $SAMPLE..."
fastqc --noextract -o "$OUTPUT_DIR" "$FASTQ_R1" "$FASTQ_R2"

echo "FastQC completed for $SAMPLE. Results saved in $OUTPUT_DIR."

