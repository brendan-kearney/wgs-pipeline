#!/bin/bash
#SBATCH --job-name=VEP-splice
#SBATCH --mem=400G
#SBATCH --output=%x-%A_%a.out

# Define paths for VEP, genome data, and sample directory
VEPDIR=/path/to/ensembl-vep
GENOMEDIR=/path/to/genomes
SAMPLE="SAMPLE_ID"
SAMPLEDIR=/path/to/genotype_$SAMPLE

# Load required modules
module load tabix

# Activate Conda environment (modify path accordingly)
source /path/to/miniconda3/bin/activate /path/to/envs/vep_env

# Define input VCF file
INPUT=$SAMPLEDIR/$SAMPLE.nodups.vcf.gz

# Run VEP with splice-detection plugins
$VEPDIR/vep --input_file $INPUT \
       --output_file $SAMPLEDIR/$SAMPLE.VEP.vcf.gz \
       --species homo_sapiens \
       --vcf \
       --fork 4 \
       --cache \
       --minimal \
       --pick \
       --force_overwrite \
       --compress_output bgzip \
       --plugin MaxEntScan,/path/to/Plugins/MaxEntScan/fordownload \
       --plugin SpliceAI,snv=$GENOMEDIR/spliceai_scores.masked.snv.hg38.vcf.gz,indel=$GENOMEDIR/spliceai_scores.masked.indel.hg38.vcf.gz \
       --plugin SpliceRegion \
       --plugin GeneSplicer,/path/to/GeneSplicer/sources/genesplicer,/path/to/GeneSplicer/human \
       --plugin SpliceVault,file=/path/to/Plugins/SpliceVault/SpliceVault_data_GRCh38.tsv.gz \
       --plugin dbscSNV,/path/to/Plugins/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
       --plugin CADD,snv=/path/to/Plugins/CADD/whole_genome_SNVs.tsv.gz,indels=/path/to/Plugins/CADD/gnomad.genomes.r4.0.indel.tsv.gz \
       --dir_cache /path/to/ensembl-vep/CACHEDIR \
       --dir_plugins /path/to/ensembl-vep/Plugins \
       --fasta /path/to/Genomes/hg38/hg38.fa

echo "Processing sample: $SAMPLE"

# Navigate to sample directory
cd $SAMPLEDIR

# Process VEP output
bgzip -d $SAMPLE.VEP.vcf.gz
grep -v '^#' $SAMPLE.VEP.vcf > $SAMPLE.noHeader.txt
tail -n +2 $SAMPLE.noHeader.txt > $SAMPLE.finalFilter.txt
bgzip $SAMPLE.VEP.vcf

# Extract specific fields
cut -f8 $SAMPLE.finalFilter.txt > $SAMPLE.onlyINFO.txt
cut -f1 $SAMPLE.finalFilter.txt > $SAMPLE.onlyHeaders.txt
cut -f2 $SAMPLE.finalFilter.txt > $SAMPLE.onlyLoc.txt

echo "Filtering splice-related annotations"

# Extract relevant SpliceAI fields
awk -F '|' '{ for (i = NF - 24; i <= NF; i++) { printf("%s", $i); if (i < NF) printf("|"); } printf("\n"); }' $SAMPLE.onlyINFO.txt > $SAMPLE.spliceAI_only.txt

# Merge extracted fields
paste $SAMPLE.onlyHeaders.txt $SAMPLE.onlyLoc.txt $SAMPLE.spliceAI_only.txt > $SAMPLE.merged_VEP.txt

# Remove intermediate files
rm $SAMPLE.onlyHeaders.txt $SAMPLE.finalFilter.txt $SAMPLE.onlyLoc.txt $SAMPLE.spliceAI_only.txt $SAMPLE.onlyINFO.txt $SAMPLE.noHeader.txt

# Remove entries with empty splice annotation
awk -F'\t' '$3 !~ /\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|/' $SAMPLE.merged_VEP.txt > $SAMPLE.filtered_VEP.txt

echo "VEP annotation completed for $SAMPLE"

### OPTIONAL: Run SQUIRLS for Splice Variant Annotation ###
# Uncomment the following lines if SQUIRLS annotation is needed.

# CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" \
#              "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" \
#              "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")
#
# chrom="${CHROMOSOMES[$SLURM_ARRAY_TASK_ID]}"
#
# java -jar -Xms2g -Xmx300g /path/to/squirls/squirls-cli.jar \
#      annotate-vcf -d /path/to/squirls -f csv -n 10000000 \
#      $SAMPLEDIR/$SAMPLE.nodups.vcf.gz $SAMPLEDIR/$SAMPLE-squirls
#
# awk -F',' '$9 != "NaN" && $9 >= 0.2' $SAMPLEDIR/$SAMPLE-squirls.csv > $SAMPLEDIR/$SAMPLE-squirls-filtered.csv
# sort -t',' -k9,9nr $SAMPLEDIR/$SAMPLE-squirls-filtered.csv > $SAMPLEDIR/$SAMPLE-squirls-sorted.csv
#
# echo "Completed SQUIRLS annotation for $SAMPLE"
# rm $SAMPLEDIR/$SAMPLE-squirls-filtered.csv
