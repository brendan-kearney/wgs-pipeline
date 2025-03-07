#!/bin/bash
#SBATCH --job-name=trimmomatic_trim
#SBATCH --mem=16G
#SBATCH --array=13-13%1

# This script runs Trimmomatic to perform QC/trimming tasks for paired-end FASTQ files.
# Currently using just ILLUMINACLIP.
# The specific trimming settings may need adjustments based on quality checks from FastQC.
# FastQC reports help identify overrepresented sequences, adapter contamination, and quality score drop-offs,
# allowing for fine-tuning of trimming parameters (e.g., additional quality filtering, leading/trailing base trimming).

# Load Trimmomatic module if required (Uncomment if using a module system)
# module load trimmomatic

# Define paths
INPUT_DIR="/path/to/processed_fastq"
WORKDIR="/path/to/analysis_directory"
SAMPLE_FILE="$WORKDIR/sample_names.txt"
OUTPUT_DIR="$WORKDIR/trimmed_fastqc"

# Read the sample name from the file based on the SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_FILE")

# Define input FASTQ files
FASTQ_R1="${INPUT_DIR}/${SAMPLE}_R1.fastq.gz"
FASTQ_R2="${INPUT_DIR}/${SAMPLE}_R2.fastq.gz"

# Define output trimmed FASTQ files
TRIMMED_R1="${OUTPUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
TRIMMED_R2="${OUTPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
UNPAIRED_R1="${OUTPUT_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
UNPAIRED_R2="${OUTPUT_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Run Trimmomatic (only ILLUMINACLIP option)
echo "Running Trimmomatic adapter trimming for $SAMPLE..."
trimmomatic PE -phred33 \
    "$FASTQ_R1" "$FASTQ_R2" \
    "$TRIMMED_R1" "$UNPAIRED_R1" \
    "$TRIMMED_R2" "$UNPAIRED_R2" \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10

echo "Trimmomatic trimming completed for $SAMPLE. Trimmed FASTQ files saved in $OUTPUT_DIR."

