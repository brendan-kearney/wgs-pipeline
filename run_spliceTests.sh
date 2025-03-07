#!/bin/bash
#SBATCH --job-name=VEP-splice
#SBATCH --mem=400G
#SBATCH --output=%x-%A_%a.out

# Define variables
VEPDIR=/data/reddylab/Reference_Data/brendan-reference/ensembl-vep
GENOMEDIR=/data/reddylab/Reference_Data/brendan-reference/bs
SAMPLE_ID=$1  # Pass sample ID as an argument
SAMPLEDIR=path_to_VCF_files/genotype_$SAMPLE_ID

module load tabix

# VEP run section: Annotates variants using the Variant Effect Predictor (VEP)
INPUT=$SAMPLEDIR/$SAMPLE_ID.hg38.vcf.gz
$VEPDIR/vep --input_file $INPUT \
       --output_file $SAMPLEDIR/$SAMPLE_ID.VEP.vcf.gz \
       --species homo_sapiens \
       --vcf \
       --fork 4 \
       --cache \
       --minimal \
       --pick \
       --force_overwrite \
       --warning_file /data/reddylab/Reference_Data/brendan-reference/warnings.txt \
       --compress_output bgzip \
       --plugin MaxEntScan,$VEPDIR/Plugins/MaxEntScan/fordownload \
       --plugin SpliceAI,snv=$GENOMEDIR/spliceai_scores.masked.snv.hg38.vcf.gz,indel=$GENOMEDIR/spliceai_scores.masked.indel.hg38.vcf.gz \
       --plugin SpliceRegion \
       --plugin GeneSplicer,$VEPDIR/GeneSplicer/sources/genesplicer,$VEPDIR/GeneSplicer/human \
       --plugin SpliceVault,file=$VEPDIR/Plugins/SpliceVault/SpliceVault_data_GRCh38.tsv.gz \
       --plugin dbscSNV,$VEPDIR/Plugins/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
       --plugin CADD,snv=$VEPDIR/Plugins/CADD/whole_genome_SNVs.tsv.gz,indels=$VEPDIR/Plugins/CADD/gnomad.genomes.r4.0.indel.tsv.gz \
       --dir_cache $VEPDIR/CACHEDIR \
       --dir_plugins $VEPDIR/Plugins \
       --fasta /data/reddylab/Reference_Data/Genomes/hg38/hg38.fa \
       --force_overwrite

# Processing VEP results
cd $SAMPLEDIR

bgzip -d $SAMPLE_ID.VEP.vcf.gz
grep -v '^#' $SAMPLE_ID.VEP.vcf > $SAMPLE_ID.noHeader.txt
tail -n +2 $SAMPLE_ID.noHeader.txt > $SAMPLE_ID.finalFilter.txt
bgzip $SAMPLE_ID.VEP.vcf

cut -f8 $SAMPLE_ID.finalFilter.txt > $SAMPLE_ID.onlyINFO.txt
cut -f1 $SAMPLE_ID.finalFilter.txt > $SAMPLE_ID.onlyHeaders.txt
cut -f2 $SAMPLE_ID.finalFilter.txt > $SAMPLE_ID.onlyLoc.txt

echo "CHECK1"

awk -F '|' '{ for (i = NF - 24; i <= NF; i++) { printf("%s", $i); if (i < NF) printf("|"); } printf("\n"); }' $SAMPLE_ID.onlyINFO.txt > $SAMPLE_ID.spliceAI_only.txt
paste $SAMPLE_ID.onlyHeaders.txt $SAMPLE_ID.onlyLoc.txt $SAMPLE_ID.spliceAI_only.txt > $SAMPLE_ID.merged_VEP.txt

rm $SAMPLE_ID.onlyHeaders.txt $SAMPLE_ID.finalFilter.txt $SAMPLE_ID.onlyLoc.txt $SAMPLE_ID.spliceAI_only.txt $SAMPLE_ID.onlyINFO.txt $SAMPLE_ID.noHeader.txt

awk -F'\t' '$3 !~ /\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|/' $SAMPLE_ID.merged_VEP.txt > $SAMPLE_ID.filtered_VEP.txt

# SQUIRLS run section: Annotates variants related to splicing defects
java -jar -Xms2g -Xmx300g /data/reddylab/btk20/software/squirls/squirls-cli-2.0.0.jar annotate-vcf \
    -d /data/reddylab/btk20/software/squirls \
    -f csv -n 10000000 \
    $SAMPLEDIR/$SAMPLE_ID.nodups.vcf.gz $SAMPLEDIR/$SAMPLE_ID-squirls

awk -F',' '$9 != "NaN" && $9 >= 0.2' $SAMPLEDIR/$SAMPLE_ID-squirls.csv > $SAMPLEDIR/$SAMPLE_ID-squirls-filtered.csv
sort -t',' -k9,9nr $SAMPLEDIR/$SAMPLE_ID-squirls-filtered.csv > $SAMPLEDIR/$SAMPLE_ID-squirls-sorted.csv
echo "completed $SAMPLE_ID squirls"

rm $SAMPLEDIR/$SAMPLE_ID-squirls-filtered.csv
