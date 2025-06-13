# =============================================================================
# Splice Variant Annotation Integration Pipeline (R Version)
#
# This script processes and integrates splicing-related variant annotations 
# from multiple sources: VEP (including MES, dbscSNV, CADD), SQUIRLS, and SpliceAI.
# It extracts relevant information, applies thresholds to flag high-impact variants,
# and merges the results into a single summary table for comparative analysis.
# The final output ranks variants by the number of high-impact calls across tools.
# =============================================================================

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# Set base directories (modify as needed)
VEP_DIR <- "/path/to/VEP_tables"
SQUIRLS_DIR <- "/path/to/squirls"
SPLICEAI_DIR <- "/path/to/splice_files"

# Function to process VEP output
process_vep_output <- function(sample_name) {
  vep_file <- file.path(VEP_DIR, paste0(sample_name, ".filtered_VEP.txt"))
  colnames <- c("chr", "loc", "info")
  VEP_df <- read_tsv(vep_file, col_names = colnames)

  # Split 'info' column into annotation fields
  info_split <- VEP_df %>%
    separate(info, into = paste0("V", 1:27), sep = "\\|", fill = "right")

  # Assign proper column names
  newcols <- c("chr", "loc", "MES_alt", "MES_diff", "MES_ref",
               "SpliceAI_DP_AG", "SpliceAI_DP_AL", "SpliceAI_DP_DG", "SpliceAI_DP_DL",
               "SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL",
               "SpliceAI_SYMBOL", "SpliceRegion", "GeneSplicer", "SpliceVault_SpliceAI_delta",
               "SpliceVault_oot_events", "SpliceVault_max_depth", "SpliceVault_pos", "SpliceVault_count",
               "SpliceVault_type", "SpliceVault_events", "dbscSNV_ada", "dbscSNV_RF",
               "CADD_phred", "CADD_raw")

  VEP_df <- bind_cols(VEP_df %>% select(-info), info_split)
  colnames(VEP_df) <- newcols

  # Convert relevant columns to numeric
  VEP_df <- VEP_df %>%
    mutate(across(c(MES_alt, MES_diff, dbscSNV_ada, dbscSNV_RF, CADD_phred, CADD_raw), as.numeric),
           MES_outcome = case_when(
             (!is.na(MES_diff) & MES_diff > 1.15 & MES_alt < 6.2) |
             (!is.na(MES_diff) & MES_diff < 0 & MES_alt > 8.5) ~ "High",
             is.na(MES_diff) | is.na(MES_alt) ~ "None",
             TRUE ~ "Low"
           )) %>%
    filter(MES_outcome != "None")

  # Extract filtered results
  MES_df <- VEP_df %>% select(chr, loc, MES_alt, MES_diff, MES_ref, MES_outcome)
  dbsc_df <- VEP_df %>%
    filter(!is.na(dbscSNV_ada) & !is.na(dbscSNV_RF)) %>%
    mutate(dbsc_outcome = if_else(dbscSNV_ada >= 0.6 | dbscSNV_RF >= 0.6, "High", "Low")) %>%
    select(chr, loc, dbscSNV_ada, dbscSNV_RF, dbsc_outcome)
  CADD_df <- VEP_df %>%
    filter(!is.na(CADD_phred)) %>%
    mutate(CADD_outcome = if_else(CADD_phred >= 10, "High", "Low")) %>%
    select(chr, loc, CADD_phred, CADD_raw, CADD_outcome)

  # Save intermediate tables
  write_csv(MES_df, file.path(VEP_DIR, paste0("MES_df_", sample_name, ".csv")))
  write_csv(dbsc_df, file.path(VEP_DIR, paste0("dbsc_df_", sample_name, ".csv")))
  write_csv(CADD_df, file.path(VEP_DIR, paste0("CADD_df_", sample_name, ".csv")))

  list(MES_df = MES_df, dbsc_df = dbsc_df, CADD_df = CADD_df)
}

# Function to merge all annotations
merge_annotations <- function(sample_name, MES_df, dbsc_df, CADD_df) {
  # Read SQUIRLS results and flag high-impact variants
  squirls_file <- file.path(SQUIRLS_DIR, paste0(sample_name, "-squirls-sorted.csv"))
  squirls_df <- read_csv(squirls_file) %>%
    distinct(chr, loc, .keep_all = TRUE) %>%
    mutate(chr = paste0("chr", chr),
           squirls_score = as.numeric(squirls_score),
           squirls_outcome = if_else(squirls_score >= 0.6, "High", "Low")) %>%
    select(chr, loc, squirls_score, squirls_outcome)

  # Read SpliceAI output and flag high scores
  spliceai_file <- file.path(SPLICEAI_DIR, paste0(sample_name, ".spliceFreq.csv"))
  SpliceAI_df <- read_tsv(spliceai_file) %>%
    mutate(DS_MAX = pmax(DS_AG, DS_AL, DS_DG, DS_DL, na.rm = TRUE),
           SpliceAI_outcome = if_else(DS_MAX >= 0.6, "High", "Low")) %>%
    select(chr, loc, DS_MAX, SpliceAI_outcome)

  # Merge all annotations using full outer joins
  comparison_df <- reduce(list(MES_df, dbsc_df, CADD_df, squirls_df, SpliceAI_df),
                          full_join, by = c("chr", "loc"))

  # Count number of 'High' and 'Low' calls per variant
  outcome_cols <- c("MES_outcome", "dbsc_outcome", "CADD_outcome", "squirls_outcome", "SpliceAI_outcome")
  comparison_df <- comparison_df %>%
    rowwise() %>%
    mutate(count_high = sum(c_across(all_of(outcome_cols)) == "High", na.rm = TRUE),
           count_low  = sum(c_across(all_of(outcome_cols)) == "Low", na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(count_high), desc(count_low))

  # Write final results to CSV
  write_csv(comparison_df, file.path(VEP_DIR, paste0(sample_name, "_final_output.csv")))
  comparison_df
}

# Run the full pipeline
sample_name <- "383582"
vep_results <- process_vep_output(sample_name)
final_result <- merge_annotations(sample_name, vep_results$MES_df, vep_results$dbsc_df, vep_results$CADD_df)
