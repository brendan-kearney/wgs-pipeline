# =============================================================================
# Structural Variant Analysis and Annotation Pipeline (R Version)
#
# This script performs structural variant (SV) analysis across multiple samples
# by integrating output from tools such as Delly, BreakDancer, Pindel, and Lumpy using
# intansv (https://github.com/YaoLab-Bioinfo/intansv)
# =============================================================================

# Load required libraries
library(intansv)
library(openxlsx)
library(ggplot2)
library(plyr)  # Ensure plyr is loaded for llply function

# Define the base directory
base_directory <- "/path/to/pipeline/SV_analysis"
ref_directory <- "/data/reddylab/Reference_Data/brendan-reference/hg38_genome/hg38_annotation.txt"

# Define the centralized output directory
global_output_dir <- file.path(base_directory, "sv_analysis_results")

# Ensure global output directory exists
if (!dir.exists(global_output_dir)) {
  dir.create(global_output_dir, recursive = TRUE)
}

# List of sample names to process
sample_names <- c("sample1", "sample2", "sample3")
# Function to process a single sample
process_sv_analysis <- function(sample_name) {
  
  # Define paths for sample
  sample_dir <- file.path(base_directory, sample_name)
  
  # Define output folder within the centralized results directory
  sample_output_dir <- file.path(global_output_dir, sample_name)
  
  # Ensure sample-specific output directory exists
  if (!dir.exists(sample_output_dir)) {
    dir.create(sample_output_dir, recursive = TRUE)
  }
  
  # Define tool output paths
  delly.file.path <- file.path(sample_dir, "delly_out", paste0(sample_name, ".delly.vcf"))
  breakdancer.file.path <- file.path(sample_dir, "breakdancer_out", paste0(sample_name, ".breakdancer.filtered.out"))
  pindel.dir.path <- file.path(sample_dir, "pindel_out")
  lumpy.file.path <- file.path(sample_dir, "lumpy_out", paste0(sample_name, ".lumpy.vcf"))
  
  # Load structural variant data (handling missing files)
  delly <- if (file.exists(delly.file.path)) readDelly(delly.file.path) else NULL
  pindel <- if (dir.exists(pindel.dir.path)) readPindel(pindel.dir.path) else NULL
  breakdancer <- if (file.exists(breakdancer.file.path)) readBreakDancer(breakdancer.file.path) else NULL
  lumpy <- if (file.exists(lumpy.file.path)) readLumpy(lumpy.file.path) else NULL
  
  # Skip if no data is loaded
  if (all(sapply(list(delly, pindel, breakdancer, lumpy), is.null))) {
    message(paste("No valid files found for sample:", sample_name, "- Skipping..."))
    return(NULL)
  }
  
  # Merge all tools' outputs
  sv_all_methods <- methodsMerge(breakdancer, delly, pindel, lumpy)
  
  # **Save each structural variant type separately**
  for (variant_type in names(sv_all_methods)) {
    if (!is.null(sv_all_methods[[variant_type]])) {
      output_file <- file.path(sample_output_dir, paste0(sample_name, ".", variant_type, ".txt"))
      write.table(sv_all_methods[[variant_type]], file = output_file, 
                  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
      message(paste("Saved:", output_file))
    }
  }
  
  # Load genome annotation file
  anno.file.path <- file.path(ref_directory, "hg38_annotation.txt")
  if (file.exists(anno.file.path)) {
    hg38_gff <- read.table(anno.file.path, header = TRUE)
  } else {
    message(paste("Annotation file missing for sample:", sample_name, "- Skipping annotation..."))
    hg38_gff <- NULL
  }
  
  # Annotate structural variants if annotation data exists
  if (!is.null(hg38_gff)) {
    
    # Remove rows with NA in "start" or "end" columns before annotation
    sv_all_methods <- llply(sv_all_methods, function(df) {
      if (!is.null(df) && "start" %in% colnames(df) && "end" %in% colnames(df)) {
        df <- na.omit(df)  # Remove NA rows
      }
      return(df)
    })
    
    # Run annotation only on cleaned data
    sv_all_methods.anno <- llply(sv_all_methods, svAnnotation, genomeAnnotation = hg38_gff)
    
    output_xlsx <- file.path(sample_output_dir, paste0(sample_name, "_sv_all_methods.anno.xlsx"))
    write.xlsx(sv_all_methods.anno, file = output_xlsx)
    
    message(paste("Saved annotation:", output_xlsx))
  }
}

# Process all samples in sample_names list
for (sample in sample_names) {
  process_sv_analysis(sample)
}

message("All samples processed successfully.")
