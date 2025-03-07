#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-13%13
#SBATCH --job-name=ali_crom
#SBATCH --output=%x-%A_%a.out

# This script runs the GATK4 PreProcessing pipeline using Cromwell to convert 
# unmapped BAM (.ubam) files into aligned and analysis-ready BAM files.
# It generates necessary input JSON and text files for Cromwell execution.

# Define directories
OUTPUT_DIR="/path/to/pipeline_gatk"
INPUT_DIR="$OUTPUT_DIR/input_uBAMs"
GENOME_DIR="$OUTPUT_DIR/hg38_genome"
INPUTS_DIR="$OUTPUT_DIR/inputs"

# Navigate to output directory
cd "$OUTPUT_DIR"

export SINGULARITY_CACHEDIR="$OUTPUT_DIR/.singularity/cache"

# Read sample names from the text file
readarray -t SAMPLES < "$OUTPUT_DIR/sample_names.txt"
FILENAME=${SAMPLES[(($SLURM_ARRAY_TASK_ID - 1))]}

# Create .ubam text input file (contains just one line with the location of the ubam file)
cd "$INPUTS_DIR"
SUFFIX=".unmapped.bam"
> "$FILENAME.bam-input.txt"
cat <<EOT >> "$FILENAME.bam-input.txt"
$INPUT_DIR/$FILENAME$SUFFIX
EOT

# Create PreProcessing JSON file for Cromwell execution
> "$FILENAME.hg38.wgs.inputs.json"
cat <<EOT >> "$FILENAME.hg38.wgs.inputs.json"
{
  "PreProcessingForVariantDiscovery_GATK4.sample_name": "$FILENAME",
  "PreProcessingForVariantDiscovery_GATK4.ref_name": "hg38",
  "PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list": "$INPUTS_DIR/$FILENAME.bam-input.txt",
  "PreProcessingForVariantDiscovery_GATK4.unmapped_bam_suffix": ".unmapped.bam",
  "PreProcessingForVariantDiscovery_GATK4.ref_dict": "$GENOME_DIR/hg38.dict",
  "PreProcessingForVariantDiscovery_GATK4.ref_fasta": "$GENOME_DIR/hg38.fa",
  "PreProcessingForVariantDiscovery_GATK4.ref_fasta_index": "$GENOME_DIR/hg38.fa.fai",
  "PreProcessingForVariantDiscovery_GATK4.ref_alt": "$GENOME_DIR/Homo_sapiens_assembly38.fasta.64.alt",
  "PreProcessingForVariantDiscovery_GATK4.ref_sa": "$GENOME_DIR/hg38.fa.sa",
  "PreProcessingForVariantDiscovery_GATK4.ref_amb": "$GENOME_DIR/hg38.fa.amb",
  "PreProcessingForVariantDiscovery_GATK4.ref_bwt": "$GENOME_DIR/hg38.fa.bwt",
  "PreProcessingForVariantDiscovery_GATK4.ref_ann": "$GENOME_DIR/hg38.fa.ann",
  "PreProcessingForVariantDiscovery_GATK4.ref_pac": "$GENOME_DIR/hg38.fa.pac",

  "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf": "$GENOME_DIR/Homo_sapiens_assembly38.dbsnp138.vcf",
  "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf_index": "$GENOME_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
  "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_VCFs": [
    "$GENOME_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "$GENOME_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz"
  ],
  "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_indices": [
    "$GENOME_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
    "$GENOME_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
  ]
}
EOT

cd "$OUTPUT_DIR"

# Run GATK4 PreProcessing pipeline using Cromwell
java -Dconfig.file=no_sql.conf -jar cromwell-78.jar run -i "$INPUTS_DIR/$FILENAME.hg38.wgs.inputs.json" "$OUTPUT_DIR/gatk4-data-processing-2.1.1/processing-for-variant-discovery-gatk4.wdl"

echo "Processing completed for $FILENAME."

