#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=geno_vcfs
#SBATCH --array=1-13%13

# This script performs joint genotyping on pre-called GVCFs using GATK4 GenotypeGVCFs.
# It then applies Variant Quality Score Recalibration (VQSR) for SNPs and indels.
# The final recalibrated VCF file is generated, along with variant statistics.
# Required tools: gatk4 and bcftools. Conda/mamba works fine.

# Define directories
RUNDIR="/path/to/pipeline_gatk"
WORKDIR="$RUNDIR/genotyped_VCFs"
GENOMEDIR="$RUNDIR/hg38_genome"
VCF_STATS_DIR="$WORKDIR/VCF_stats"

# Read sample names from the file
readarray -t SAMPLES < "$RUNDIR/sample_names.txt"
COHORTNAME=${SAMPLES[(($SLURM_ARRAY_TASK_ID - 1))]}

# Ensure necessary directories exist
mkdir -p "$WORKDIR/genotype_$COHORTNAME"
mkdir -p "$VCF_STATS_DIR"

cd "$WORKDIR"

# Perform joint genotyping using GATK GenotypeGVCFs
echo "Running joint genotyping for $COHORTNAME..."
gatk GenotypeGVCFs \
	-R "$GENOMEDIR/hg38.fa" \
	-V "$RUNDIR/finished_GVCFs/$COHORTNAME.hg38.g.vcf.gz" \
	-O "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.genotyped.vcf.gz"

# Recalibrate indels
echo "Performing VQSR for indels on $COHORTNAME..."
gatk VariantRecalibrator \
	-V "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.genotyped.vcf.gz" \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
	-mode INDEL \
	--max-gaussians 4 \
	-resource:mills,known=false,training=true,truth=true,prior=12 "$GENOMEDIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
	-resource:axiomPoly,known=false,training=true,truth=false,prior=10 "$GENOMEDIR/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz" \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2 "$GENOMEDIR/Homo_sapiens_assembly38.dbsnp138.vcf" \
	-O "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.indels.recal" \
	--tranches-file "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.indels.tranches"

# Recalibrate SNPs
echo "Performing VQSR for SNPs on $COHORTNAME..."
gatk VariantRecalibrator \
	-V "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.genotyped.vcf.gz" \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
	-mode SNP \
	--max-gaussians 6 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15 "$GENOMEDIR/hapmap_3.3.hg38.vcf.gz" \
	-resource:omni,known=false,training=true,truth=true,prior=12 "$GENOMEDIR/1000G_omni2.5.hg38.vcf.gz" \
	-resource:1000G,known=false,training=true,truth=false,prior=10 "$GENOMEDIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
	-resource:dbsnp,known=true,training=false,truth=false,prior=7 "$GENOMEDIR/Homo_sapiens_assembly38.dbsnp138.vcf" \
	-O "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.snp.recal" \
	--tranches-file "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.snp.tranches"

# Apply VQSR for indels
echo "Applying VQSR for indels on $COHORTNAME..."
gatk ApplyVQSR \
	-V "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.genotyped.vcf.gz" \
	--recal-file "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.indels.recal" \
	--tranches-file "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.indels.tranches" \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	--mode INDEL \
	-O "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.indel.recalibrated.vcf.gz"

# Apply VQSR for SNPs on the indel-filtered VCF
echo "Applying VQSR for SNPs on $COHORTNAME..."
gatk ApplyVQSR \
	-V "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.indel.recalibrated.vcf.gz" \
	--recal-file "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.snp.recal" \
	--tranches-file "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.snp.tranches" \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	--mode SNP \
	-O "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.recalibrated.vcf.gz"

# Generate statistics using bcftools
echo "Generating VCF statistics for $COHORTNAME..."
bcftools stats "$WORKDIR/genotype_$COHORTNAME/$COHORTNAME.recalibrated.vcf.gz" > "$VCF_STATS_DIR/$COHORTNAME.vcf.stats"

# Cleanup intermediate files
echo "Cleaning up intermediate files for $COHORTNAME..."
cd "$WORKDIR/genotype_$COHORTNAME"
rm "$COHORTNAME.indels.recal"
rm "$COHORTNAME.indels.recal.idx"
rm "$COHORTNAME.indels.tranches"
rm "$COHORTNAME.snp.recal"
rm "$COHORTNAME.snp.recal.idx"
rm "$COHORTNAME.snp.tranches"
rm "$COHORTNAME.indel.recalibrated.vcf.gz"
rm "$COHORTNAME.indel.recalibrated.vcf.gz.tbi"
rm "$COHORTNAME.genotyped.vcf.gz"
rm "$COHORTNAME.genotyped.vcf.gz.tbi"

echo "Genotyping and VQSR completed for $COHORTNAME."

