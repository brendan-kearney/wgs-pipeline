#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --array=13-13%1
#SBATCH --job-name=fastq_to_ubam
#SBATCH --output=%x_%A_%a.out

# This script converts trimmed FASTQ files into unmapped BAM (.u.bam) format for alignment preparation.
# It builds an input JSON configuration for the GATK seq-format-conversion workflow and runs Cromwell.
# This step is essential for maintaining read group information and proper alignment downstream.
# The extra json info (library name, platform, etc) are unnecessary to run the script

# Define directories
SAMPLEDIR="/path/to/trimmed_fastqc"
WORKDIR="/path/to/pipeline_gatk"
CROMWELLDIR="/path/to/cromwell_jar_file"

cd "$WORKDIR"

export SINGULARITY_CACHEDIR="$WORKDIR/.singularity/cache"

# Read sample names from a text file
readarray -t SAMPLES < "$WORKDIR/sample_names.txt"
FILENAME=${SAMPLES[(($SLURM_ARRAY_TASK_ID - 1))]}

mkdir -p inputs
cd inputs

R1_SUFFIX="_R1_trimmed.fastq.gz"
R2_SUFFIX="_R2_trimmed.fastq.gz"

# Build JSON input file for each sample
> "$FILENAME.format-conversion.json"
cat <<EOT >> "$FILENAME.format-conversion.json"
{
  "ConvertPairedFastQsToUnmappedBamWf.readgroup_name": "$FILENAME",
  "ConvertPairedFastQsToUnmappedBamWf.sample_name": "$FILENAME",
  "ConvertPairedFastQsToUnmappedBamWf.fastq_1": "$SAMPLEDIR/$FILENAME$R1_SUFFIX",
  "ConvertPairedFastQsToUnmappedBamWf.fastq_2": "$SAMPLEDIR/$FILENAME$R2_SUFFIX",
  "ConvertPairedFastQsToUnmappedBamWf.library_name": "default_library",
  "ConvertPairedFastQsToUnmappedBamWf.platform_unit": "default_unit",
  "ConvertPairedFastQsToUnmappedBamWf.run_date": "2025-01-01T00:00:00+0000",
  "ConvertPairedFastQsToUnmappedBamWf.platform_name": "illumina",
  "ConvertPairedFastQsToUnmappedBamWf.sequencing_center": "default_center",
  "ConvertPairedFastQsToUnmappedBamWf.make_fofn": true  
}
EOT

# Execute WDL script with Cromwell
cd "$WORKDIR"
java -Dconfig.file=no_sql.conf -jar cromwell-78.jar run -i "$WORKDIR/inputs/$FILENAME.format-conversion.json" seq-format-conversion/paired-fastq-to-unmapped-bam.wdl

echo "Processing completed for $FILENAME."

