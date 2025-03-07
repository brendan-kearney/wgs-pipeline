# Whole Genome Sequencing Variant Calling Pipeline

## Overview

This pipeline processes whole genome sequencing (WGS) data from raw FASTQ files to genotyped and recalibrated variant call format (VCF) files. The workflow follows these main steps:

1. Concatenation of FASTQ files from multiple sequencing lanes.
2. Quality control (QC) checks using FastQC.
3. Adapter trimming with Trimmomatic.
4. Conversion of FASTQ to uBAM format for alignment.
5. Preprocessing of reads, including alignment and duplicate marking.
6. Variant calling to produce GVCFs.
7. Genotyping of individual samples.
8. Joint genotyping for related individuals (families).

Each step utilizes specific tools and reference files, ensuring high-quality variant detection.

---

## Required Tools & Dependencies

| Tool       | Version  | Purpose                                      | Used in Script(s) |
|------------|---------|----------------------------------------------|-------------------|
| FastQC     | Any recent | Quality control of FASTQ files               | `2-run_fastqc.sh` |
| Trimmomatic| Any recent | Adapter removal and quality trimming         | `3-run_trimmomatic.sh` |
| Cromwell   | 78+     | WDL Workflow Engine                          | `4-cromwell_seqFormat.sh`, `5-cromwell_preProcessing.sh`, `6-cromwell_variantCalling.sh` |
| GATK       | 4.1.0.0 | Preprocessing, Variant Calling, and Genotyping | `5-cromwell_preProcessing.sh`, `6-cromwell_variantCalling.sh`, `7-genotypeSingleSamples.sh`, `8-genotypeFamilies.sh` |
| Samtools   | Any recent | BAM file indexing & statistics               | `6-cromwell_variantCalling.sh` |
| BCFTools   | Any recent | Variant statistics for VCFs                   | `8-genotypeFamilies.sh` |
| Python/R   | Optional | Additional post-processing & statistics       | Not explicitly mentioned |

---

## Reference Files Required

Most hg38 references can be found online via 1000G and GATK resources: [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

The `hg38_wgs_custom_intervals.txt` file contains the path to `wgs_calling_regions.hg38.interval_list`.

A complete hg38 reference directory is available at:
```
/data/reddylab/Reference_Data/brendan-reference/hg38_genome
```

---

## Pipeline Step Descriptions

### Step 1: Concatenation of Multi-Lane FASTQ Files
- **Script**: `1-concat_fastqs.sh`
- **Tools**: Bash (no additional dependencies)
- **Purpose**: Merges multiple sequencing lanes into a single FASTQ file per read strand (`_R1.fastq.gz` and `_R2.fastq.gz`).
- **Input**: Raw sequencing lane FASTQ files (`Sample_L001_R1_001.fastq.gz`, etc.)
- **Output**: Merged FASTQ files (`sampleA_R1.fastq.gz`, `sampleA_R2.fastq.gz`)

### Step 2: FastQC Quality Control
- **Script**: `2-run_fastqc.sh`
- **Tools**: FastQC
- **Purpose**: Performs quality control checks on raw FASTQ files.
- **Output**: FastQC reports (`sampleA_fastqc.html`, `.zip`)

### Step 3: Adapter Trimming with Trimmomatic
- **Script**: `3-run_trimmomatic.sh`
- **Tools**: Trimmomatic
- **Purpose**: Removes Illumina adapters and generates trimmed paired-end FASTQ files.
- **Output**: Trimmed FASTQ files (`sampleA_R1_trimmed.fastq.gz`, `sampleA_R2_trimmed.fastq.gz`)

### Step 4: Convert FASTQ to Unmapped BAM (uBAM)
- **Script**: `4-cromwell_seqFormat.sh`
- **Tools**: Cromwell
- **Purpose**: Converts FASTQ to uBAM format, preparing it for alignment.
- **Output**: Unmapped BAM file (`sampleA.unmapped.bam`)

### Step 5: Read Preprocessing & Alignment
- **Script**: `5-cromwell_preProcessing.sh`
- **Tools**: Cromwell, GATK
- **Purpose**: Aligns uBAM files to the reference genome (`hg38.fa`), marks duplicate reads, and sorts BAM files.
- **Output**: Aligned BAM files (`sampleA.hg38.bam`)

### Step 6: Variant Calling with HaplotypeCaller
- **Script**: `6-cromwell_variantCalling.sh`
- **Tools**: Cromwell, GATK, Samtools
- **Purpose**: Calls variants from aligned BAM files using HaplotypeCaller and generates GVCF files per sample.
- **Output**: GVCF file (`sampleA.hg38.g.vcf.gz`)

### Step 7: Genotyping Individual Samples
- **Script**: `7-genotypeSingleSamples.sh`
- **Tools**: GATK
- **Purpose**: Converts GVCF files into genotyped VCF files and applies Variant Quality Score Recalibration (VQSR).
- **Output**: Recalibrated VCF (`sampleA.recalibrated.vcf.gz`)

### Step 8: Joint Genotyping for Families
- **Script**: `8-genotypeFamilies.sh`
- **Tools**: GATK, BCFTools
- **Purpose**: Performs joint genotyping across multiple related individuals and applies VQSR.
- **Output**:
  - Family-level recalibrated VCF (`sample_family.recalibrated.vcf.gz`)
  - Variant statistics (`sample_family.vcf.stats`)

---

## How to Run

The pipeline uses `no_sql.conf` as a Cromwell configuration file, which enables a workaround for Docker on HPC clusters where Docker is disabled. Ensure that `CACHE_DIR` is properly set (lines 74 and 76).

Cromwell generates outputs in `cromwell_executions`, structured as:
```
cromwell_executions/HaplotypeCallerGvcf_GATK4/<unique_id>/call-HaplotypeCaller/shard-0/execution/sample.hg38.g.vcf
```
These files will need manual processing before the next step.

### Estimated Runtime
- **Whole Genome Sequencing (WGS)**: ~1-5 days (depends on processing power and sample size).
- **Whole Exome Sequencing (WES)**: Much shorter.

---

---

## EXOMISER - How to Run

Exomiser is a Java program that identifies potential disease-causing variants from whole-exome or whole-genome sequencing data using phenotype-driven algorithms.

### Steps of an Exomiser Run:
1. **Inputs**: A joint-called VCF trio file and an associated `.ped` file.
2. **Annotation** of variants using Jannovar.
3. **Removal** of variants based on user-supplied criteria.
4. **Phenotypic assessment** and ranking of variants by predicted pathogenicity.

### Filtering Options:
- **Failed previous VCF filters**: Removes variants failing GATK HaplotypeCaller and VQSR filtering.
- **Regulatory feature**: Removes non-regulatory, non-coding variants >20Kb from a known gene.
- **Frequency**: Removes variants with frequency >2.0% in any source (gnomAD, ExAC, ESP, DBSNP).
- **Pathogenicity**: Removes variants with Exomiser Variant Score < 0.500.
- **Inheritance-based MAF thresholds**:
  - **Autosomal Dominant**: 0.1%
  - **Autosomal Recessive (Homozygous Alt)**: 0.1%
  - **Autosomal Recessive (Compound Heterozygous)**: 2.0%
  - **X-Linked Dominant**: 0.1%
  - **X-Linked Recessive (Homozygous Alt)**: 0.1%
  - **X-Linked Recessive (Compound Heterozygous)**: 2.0%
- **Phenotypes**: Filters for genes associated with provided HPO terms.

### Exomiser Output Scores:
1. **Variant Score** - Based on rarity and pathogenicity predictions.
2. **Gene Variant Score** - Averaged across passed variants of a gene.
3. **Gene-Specific Phenotype Score** - Matches patient HPO terms to known disease annotations.
4. **Exomiser Score** - Combined phenotype and gene variant score.

All Exomiser analyses should be run in the `EXOMISER` directory, requiring `.yaml` and `.ped` files, `application.properties`, and `exomiser.jar` (version 14.0.0).

---

## Splice and Structural Variant Detection

Both operate with two scripts, starting with a script to run all testing tools in parallel (run_spliceTests.sh and sv_analysis.sh), and an integration script to compile data and output results (splice_tools.sh and SV_processing.R). 
Structural variant analysis is in its own folder.

---

## Acknowledgments

This pipeline leverages several open-source tools and repositories:

- **[GATK workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)** – Base pipeline for short variant discovery.
- **[GATK](https://github.com/broadinstitute/gatk)** – For variant calling, preprocessing, and genotyping.
- **[Cromwell](https://github.com/broadinstitute/cromwell)** – Workflow execution engine for running WDL pipelines.
- **[Samtools](https://github.com/samtools/samtools)** – BAM file manipulation and indexing.
- **[BCFTools](https://github.com/samtools/bcftools)** – Variant filtering and statistics.
- **[FastQC](https://github.com/s-andrews/FastQC)** – Quality control of sequencing data.
- **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** – Adapter and quality trimming for FASTQ reads.
- **[Exomiser](https://github.com/exomiser/Exomiser)** – Phenotype-driven variant prioritization.
- **[ensembl-vep](https://github.com/Ensembl/ensembl-vep)** – Genomic variant annotator.
- **[SQUIRLS](https://github.com/monarch-initiative/Squirls)** – Splice prediction.
- **[breakdancer](https://github.com/genome/breakdancer)** - Structural Variant detection.
- **[LUMPY](https://github.com/arq5x/lumpy-sv)** – Structural variant detection.
- **[Delly](https://github.com/dellytools/delly)** – Structural variant discovery.
- **[Pindel](https://github.com/genome/pindel)** – Indel detection in NGS data.
- **[CNVNator](https://github.com/abyzovlab/CNVnator)** - CNV discovery.
- **[intansv](https://github.com/YaoLab-Bioinfo/intansv)** - Integrative analysis of Structural Variants
