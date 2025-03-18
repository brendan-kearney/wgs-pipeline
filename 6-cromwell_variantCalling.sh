#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-13%13
#SBATCH --job-name=crom_var
#SBATCH --output=%x-%A_%a.out

# This script runs variant calling using GATK4 HaplotypeCaller in GVCF mode via Cromwell.
# It performs BAM indexing, computes BAM statistics using samtools, 
# and generates a JSON input file for HaplotypeCaller execution
# BAM stat file is useful in further analysis (structural variants)

module load samtools

# Define directories
WORKDIR="/path/to/pipeline_gatk"
BAM_DIR="$WORKDIR/finished_bams"
STATS_DIR="$BAM_DIR/bam_stats"
GENOME_DIR="$WORKDIR/hg38_genome"
INPUTS_DIR="$WORKDIR/inputs"
CROMWELLDIR="/path/to/cromwell_jar_file"

# Navigate to working directory
cd "$WORKDIR"

export SINGULARITY_CACHEDIR="$WORKDIR/.singularity/cache"

# Read sample names from the text file
readarray -t SAMPLES < "$WORKDIR/sample_names.txt"
FILENAME=${SAMPLES[(($SLURM_ARRAY_TASK_ID - 1))]}

echo "Processing sample: $FILENAME"

# Ensure necessary directories exist
mkdir -p "$STATS_DIR"

# Compute BAM statistics
echo "Calculating BAM statistics for $FILENAME..."
samtools stats "$BAM_DIR/$FILENAME.hg38.bam" > "$STATS_DIR/$FILENAME.bam.stats"

# Update the intervals file dynamically
echo "Updating interval file..."
cat "$GENOME_DIR/wgs_calling_regions.hg38.interval_list" | tr '\n' ' ' > "$GENOME_DIR/hg38_wgs_custom_intervals.txt"
echo "" >> "$GENOME_DIR/hg38_wgs_custom_intervals.txt"  # Ensure the file ends with a newline

# Generate JSON input file for Cromwell execution
echo "Generating JSON input for HaplotypeCaller..."
cd "$INPUTS_DIR"
> "$FILENAME.haplotypecaller-inputs.json"
cat <<EOT >> "$FILENAME.haplotypecaller-inputs.json"
{
  "HaplotypeCallerGvcf_GATK4.input_bam": "$BAM_DIR/$FILENAME.hg38.bam",
  "HaplotypeCallerGvcf_GATK4.input_bam_index": "$BAM_DIR/$FILENAME.hg38.bai",
  "HaplotypeCallerGvcf_GATK4.ref_dict": "$GENOME_DIR/hg38.dict",
  "HaplotypeCallerGvcf_GATK4.ref_fasta": "$GENOME_DIR/hg38.fa",
  "HaplotypeCallerGvcf_GATK4.ref_fasta_index": "$GENOME_DIR/hg38.fa.fai",
  "HaplotypeCallerGvcf_GATK4.scattered_calling_intervals_list": "$GENOME_DIR/hg38_wgs_custom_intervals.txt"
}
EOT

# Run Cromwell with the generated JSON input
echo "Starting Cromwell for $FILENAME..."
cd "$WORKDIR"
java -Dconfig.file=no_sql.conf -jar $CROMWELLDIR/cromwell-78.jar run -i "$INPUTS_DIR/$FILENAME.haplotypecaller-inputs.json" gatk4-germline-snps-indels/haplotypecaller-gvcf-gatk4.wdl

echo "Variant calling completed for $FILENAME."

