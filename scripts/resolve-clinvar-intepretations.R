################################################################################
# resolve-clinvar-intepretations.R
# written by Ryan Corbett, Patricia Sullivan
#
# This script selects unique ClinVar variant submission calls based on a list of
# predetermined criteria, to be used in AutoGVP for ClinVar variants that need
# resolving due to conflicting calls
#
# usage: resolve-clinvar-intepretations.R --variant_summary <variant file>
#                                         --submission_summary <submission file>
#
# NOTE: this script must be run BEFORE running run_autogvp.sh
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
  library("lubridate")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# parse parameters
option_list <- list(
  make_option(c("--variant_summary"),
    type = "character",
    help = "variant_summary file (format: variant_summary_YYYY-MM.txt.gz)"
  ),
  make_option(c("--submission_summary"),
    type = "character",
    help = "specific submission summary file (format: submission_summary_YYYY-MM.txt.gz)"
  ),
  make_option(c("--outdir"),
    type = "character",
    help = "output directory"
  ),
  make_option(c("--conceptID_list"),
    type = "character",
    default = NULL,
    help = "list of ClinVar concept IDs"
  ),
  make_option(c("--conflict_res"),
    type = "character",
    default = "latest",
    help = "how to resolve conflicting interpretations? 'latest' or 'most_severe'"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (read)
input_submission_file <- opt$submission_summary
input_variant_summary <- opt$variant_summary
results_dir <- opt$outdir
conceptID_file <- opt$conceptID_list
conflict_res <- opt$conflict_res

# Function to extract date (YYYY-MM) from filename
extract_date <- function(path) {
  parts <- strsplit(basename(path), "_")[[1]]
  last_part <- tail(parts, 1)
  sub("\\.tsv\\.gz$|\\.txt\\.gz$", "", last_part)
}

# Extract dates
date_submission <- extract_date(input_submission_file)
date_variant <- extract_date(input_variant_summary)

# Check if they match
if (date_submission != date_variant) {
  stop(sprintf("ClinVar date mismatch: %s vs %s", date_submission, date_variant))
}

# Assign concept label to be overwritten if used
concept_label <- ""

## load variant summary file, which reports latest ClinVar consensus calls for each variant
variant_summary_df <- vroom(input_variant_summary,
  delim = "\t",
  col_types = c(ReferenceAlleleVCF = "c", AlternateAlleleVCF = "c", PositionVCF = "i", VariationID = "n"),
  show_col_types = FALSE
) %>%
  dplyr::rename(
    AlleleID = dplyr::any_of(c("AlleleID", "#AlleleID"))
  ) %>%
  dplyr::filter(
    Assembly == "GRCh38", # retain only variants mapped to hg38
    ReferenceAlleleVCF != "na",
    AlternateAlleleVCF != "na"
  ) %>%
  # add vcf id column
  dplyr::mutate(
    vcf_id = str_remove_all(paste(Chromosome, "-", PositionVCF, "-", ReferenceAlleleVCF, "-", AlternateAlleleVCF), " "),
    vcf_id = str_replace_all(vcf_id, "chr", ""),
    VariationID = as.double(noquote(VariationID)),
    LastEvaluated = case_when(
      LastEvaluated == "-" ~ NA_character_,
      TRUE ~ LastEvaluated
    )
  ) %>%
  dplyr::filter(!ReviewStatus %in% c(
    "no assertion provided",
    "no assertion criteria provided",
    "no classification for the individual variant",
    "no classification provided",
    "no classification for the single variant",
    "no classifications from unflagged records"
  )) %>%
  dplyr::mutate(
    ClinSig_resolved = case_when(
      grepl("Conflicting classifications of pathogenicity", ClinicalSignificance) ~ "Conflicting classifications of pathogenicity",
      grepl("Uncertain significance", ClinicalSignificance) ~ "Uncertain significance",
      grepl("Pathogenic/Likely pathogenic", ClinicalSignificance) ~ "Pathogenic/Likely pathogenic",
      grepl("Benign/Likely benign", ClinicalSignificance) ~ "Benign/Likely benign",
      grepl("Likely pathogenic", ClinicalSignificance) ~ "Likely pathogenic",
      grepl("Likely benign", ClinicalSignificance) ~ "Likely benign",
      grepl("Pathogenic", ClinicalSignificance) ~ "Pathogenic",
      grepl("Benign", ClinicalSignificance) ~ "Benign",
      grepl("Uncertain risk allele", ClinicalSignificance) ~ "Uncertain significance",
      grepl("risk allele", ClinicalSignificance) ~ "Risk allele",
      grepl("VUS", ClinicalSignificance) ~ "Uncertain significance",
      TRUE ~ NA
    ),
    Stars = case_when(
      ReviewStatus == "practice guideline" ~ 4,
      ReviewStatus == "reviewed by expert panel" ~ 3,
      ReviewStatus == "criteria provided, multiple submitters, no conflicts" ~ 2,
      ReviewStatus %in% c("criteria provided, conflicting classifications", "criteria provided, single submitter") ~ 1,
      TRUE ~ NA
    )
  ) %>%
  dplyr::filter(!is.na(ClinSig_resolved))


# Load ClinVar submission summary file, which reports all submissions for each ClinVar variant

# Open the file and read lines until a line without '#' is found
con <- file(input_submission_file, "r")
skip_lines <- -1 # keep last line for column names
while (grepl("^#", readLines(con, n = 1))) {
  skip_lines <- skip_lines + 1
}
close(con)
skip_lines <- max(skip_lines, 0)

# Load submission file while skipping number of lines determined above
submission_summary_df <- vroom(input_submission_file,
  skip = skip_lines,
  delim = "\t",
  quote = "",
  show_col_types = F
) %>%
  dplyr::rename(
    VariationID = dplyr::any_of(c("VariationID", "#VariationID"))
  )

# Define columns to check the presence of
required_cols <- c(
  "VariationID",
  "ClinicalSignificance",
  "DateLastEvaluated",
  "ReviewStatus"
)

# Columns present in submissions file as of 01 2026
#  "#VariationID", "ClinicalSignificance", "DateLastEvaluated",
#  "Description", "SubmittedPhenotypeInfo", "ReportedPhenotypeInfo",
#  "ReviewStatus", "CollectionMethod", "OriginCounts", "Submitter",
#  "SCV", "SubmittedGeneSymbol", "ExplanationOfInterpretation",
#  "SomaticClinicalImpact",	"Oncogenicity",	"ContributesToAggregateClassification"

missing_cols <- setdiff(required_cols, names(submission_summary_df))
if (length(missing_cols) > 0) {
  stop(
    paste(
      "Required column(s) missing from ClinVar submission summary file:",
      paste(missing_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (!"ContributesToAggregateClassification" %in% names(submission_summary_df)) {
  message(
    "Note: 'ContributesToAggregateClassification' column not found.",
    "This is likely an older ClinVar submission file.",
    "Proceeding without aggregate classification filtering."
  )
}

submission_summary_df <- submission_summary_df %>%
  # Redefine `DateLastEvaluated`; remove all quote characters from Description
  dplyr::mutate(
    DateLastEvaluated = case_when(
      DateLastEvaluated == "-" ~ NA_character_,
      TRUE ~ DateLastEvaluated
    ),
    Description = str_remove_all(Description, "[\"'‘’“”]"),
    VariationID = as.double(VariationID)
  ) %>%
  dplyr::filter(
    !ReviewStatus %in% c(
      "no assertion provided",
      "no assertion criteria provided",
      "no classification provided",
      "flagged submission"
    ),
    ClinicalSignificance %in% c(
      "Pathogenic",
      "Likely pathogenic",
      "Pathogenic, low penetrance",
      "Likely pathogenic, low penetrance",
      "Established risk allele",
      "Likely risk allele",
      "Benign",
      "Likely benign",
      "Uncertain significance",
      "Uncertain risk allele",
      "VUS-high",
      "VUS-mid",
      "VUS-low"
    ),
    # Filter on contributing records if this column is present in supplied file
    if ("ContributesToAggregateClassification" %in% names(.)) {
      ContributesToAggregateClassification == "yes"
    } else {
      TRUE
    }
  )

# merge submission_summary and variant_summary info
submission_merged_df <- submission_summary_df %>%
  dplyr::rename("LastEvaluated" = DateLastEvaluated) %>%
  inner_join(variant_summary_df,
    by = "VariationID",
    multiple = "all",
    suffix = c("_sub", "_var"),
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    LastEvaluated = coalesce(LastEvaluated_sub, LastEvaluated_var),
    ClinSig_consensus = case_when(
      ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic", "Pathogenic, low penetrance", "Likely pathogenic, low penetrance", "Established risk allele", "Likely risk allele") ~ "Pathogenic/Likely pathogenic",
      ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "Benign/Likely benign",
      ClinicalSignificance_sub %in% c("Uncertain significance", "Uncertain risk allele", "VUS-high", "VUS-mid", "VUS-low") ~ "Uncertain significance",
      TRUE ~ NA_character_
    )
  )

# Extract submissions that have no conflicts
variants_resolved <- submission_merged_df %>%
  filter(
    ClinSig_resolved != "Conflicting classifications of pathogenicity",
    !is.na(ClinSig_resolved)
  ) %>%
  dplyr::arrange(desc(mdy(LastEvaluated))) %>%
  distinct(VariationID, .keep_all = T) %>%
  mutate(ClinSig_report = ClinSig_resolved) %>%
  arrange(VariationID)

# IF list of concept IDs provided -- filter remaining submissions to only those associated with concept IDs, and resolve conflicts by consensus, latest date, or severity
if (!is.null(conceptID_file)) {
  # pull out concept ID file category (cancer, all, cpg)
  concept_cat <- strsplit(basename(conceptID_file), "_")[[1]][2]
  concept_label <- glue::glue("-{concept_cat}-{conflict_res}")

  # read in concept IDs
  conceptIDs <- read_lines(conceptID_file)

  # Filter variants for those associated with concept IDs
  variants_with_conceptIDs <- submission_merged_df %>%
    dplyr::filter(ClinSig_resolved == "Conflicting classifications of pathogenicity") %>%
    dplyr::mutate(conceptID = unlist(lapply(strsplit(ReportedPhenotypeInfo, ":"), function(x) x[[1]]))) %>%
    dplyr::filter(conceptID %in% conceptIDs) %>%
    dplyr::select(-conceptID)

  # if no. submissions remaining = 1, add to no conflict variants
  variants_no_conflicts_conceptID <- variants_with_conceptIDs %>%
    count(VariationID) %>%
    dplyr::filter(n == 1) %>%
    pull(VariationID)

  # add variants with one submission to variants_no_conflicts
  variants_resolved <- variants_with_conceptIDs %>%
    filter(VariationID %in% variants_no_conflicts_conceptID) %>%
    dplyr::mutate(
      ClinSig_report = ClinicalSignificance_sub,
      ReviewStatus_var = glue::glue("{ReviewStatus_var}, single submission associated with conceptIDs")
    ) %>%
    bind_rows(variants_resolved)

  # Compute max-count buckets per conceptID variant once
  clinsig_counts_conceptIDs <- variants_with_conceptIDs %>%
    filter(!VariationID %in% variants_no_conflicts_conceptID) %>%
    dplyr::filter(!is.na(ClinSig_consensus)) %>%
    count(VariationID, ClinSig_consensus) %>%
    group_by(VariationID) %>%
    dplyr::filter(n == max(n))

  consensus_calls_conceptIDs <- clinsig_counts_conceptIDs %>%
    dplyr::filter(n() == 1) %>%
    ungroup()
  tied_calls_conceptIDs <- clinsig_counts_conceptIDs %>%
    dplyr::filter(n() > 1) %>%
    ungroup()

  # extract variants resolved through consensus calling
  variants_resolved <- variants_with_conceptIDs %>%
    dplyr::semi_join(consensus_calls_conceptIDs, by = c("VariationID", "ClinSig_consensus")) %>%
    dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
    distinct(VariationID, .keep_all = T) %>%
    arrange(VariationID) %>%
    dplyr::mutate(ReviewStatus_var = glue::glue("{ReviewStatus_var}, submissions associated with conceptIDs, consensus call taken")) %>%
    dplyr::rename(ClinSig_report = ClinSig_consensus) %>%
    bind_rows(variants_resolved)

  # Resolve remaining conflicts for variants associated with concept IDs, either by taking date last evaluated or most severe call
  if (conflict_res == "latest") {
    variants_resolved <- variants_with_conceptIDs %>%
      dplyr::filter(
        !VariationID %in% consensus_calls_conceptIDs$VariationID,
        !VariationID %in% variants_no_conflicts_conceptID
      ) %>%
      dplyr::semi_join(tied_calls_conceptIDs, by = c("VariationID", "ClinSig_consensus")) %>%
      dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
      distinct(VariationID, .keep_all = T) %>%
      dplyr::mutate(
        ClinSig_report = ClinicalSignificance_sub,
        ReviewStatus_var = glue::glue("{ReviewStatus_var}, submissions associated with conceptIDs, latest date evaluated taken")
      ) %>%
      bind_rows(variants_resolved)
  }

  if (conflict_res == "most_severe") {
    variants_resolved <- variants_with_conceptIDs %>%
      dplyr::filter(
        !VariationID %in% consensus_calls_conceptIDs$VariationID,
        !VariationID %in% variants_no_conflicts_conceptID
      ) %>%
      dplyr::semi_join(tied_calls_conceptIDs, by = c("VariationID", "ClinSig_consensus")) %>%
      dplyr::mutate(ClinicalSignificance_sub = fct_relevel(
        ClinicalSignificance_sub,
        c(
          "Pathogenic", "Likely pathogenic",
          "Pathogenic, low penetrance", "Likely pathogenic, low penetrance",
          "Established risk allele", "Likely risk allele",
          "VUS-high", "VUS-mid", "Uncertain significance", "Uncertain risk allele", "VUS-low",
          "Likely benign", "Benign"
        )
      )) %>%
      dplyr::arrange(VariationID, ClinicalSignificance_sub, desc(mdy(LastEvaluated_sub))) %>%
      distinct(VariationID, .keep_all = T) %>%
      dplyr::mutate(
        ClinSig_report = ClinicalSignificance_sub,
        ReviewStatus_var = glue::glue("{ReviewStatus_var}, submissions associated with conceptIDs, most severe call at latest date taken")
      ) %>%
      bind_rows(variants_resolved)
  }
}

# Extract remaining variants with conflicts
conflicting_variants <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% variants_resolved$VariationID)

# Compute max-count buckets per variant once
clinsig_counts <- conflicting_variants %>%
  dplyr::filter(!is.na(ClinSig_consensus)) %>%
  count(VariationID, ClinSig_consensus) %>%
  group_by(VariationID) %>%
  dplyr::filter(n == max(n))

consensus_calls <- clinsig_counts %>%
  dplyr::filter(n() == 1) %>%
  ungroup()
tied_calls <- clinsig_counts %>%
  dplyr::filter(n() > 1) %>%
  ungroup()

# Extract variants with consensus calls
variants_consensus_call <- conflicting_variants %>%
  dplyr::semi_join(consensus_calls, by = c("VariationID", "ClinSig_consensus")) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  distinct(VariationID, .keep_all = T) %>%
  arrange(VariationID) %>%
  dplyr::mutate(ReviewStatus_var = glue::glue("{ReviewStatus_var}, consensus call taken")) %>%
  dplyr::rename(ClinSig_report = ClinSig_consensus)

# Resolve conflicts in remaining variants by taking call at date last evaluated
variants_conflicts_latest <- conflicting_variants %>%
  dplyr::filter(!VariationID %in% variants_consensus_call$VariationID) %>%
  dplyr::semi_join(tied_calls, by = c("VariationID", "ClinSig_consensus")) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::mutate(
    ClinSig_report = ClinicalSignificance_sub,
    ReviewStatus_var = glue::glue("{ReviewStatus_var}, latest date evaluated call taken")
  )

# create final df of all resolved submissions, and take Clinical Significance of submission as ClinSig
submission_final_df <- variants_resolved %>%
  bind_rows(variants_consensus_call, variants_conflicts_latest) %>%
  dplyr::mutate(
    ClinicalSignificance = case_when(
      grepl("Uncertain risk allele", ClinSig_report, ignore.case = TRUE) ~ "Uncertain significance",
      grepl("risk allele", ClinSig_report, ignore.case = TRUE) ~ "Risk allele",
      TRUE ~ ClinSig_report
    ),
    ReviewStatus = ReviewStatus_var,
    SubmittedPhenotypeInfo = case_when(
      SubmittedPhenotypeInfo == "Not Provided" ~ NA_character_,
      TRUE ~ SubmittedPhenotypeInfo
    ),
    Description = case_when(
      Description == "-" ~ NA_character_,
      TRUE ~ Description
    )
  ) %>%
  filter(
    !is.na(ClinicalSignificance),
    !ClinicalSignificance %in% c("-", "not provided", "association", "risk factor", "drug response")
  ) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  distinct(VariationID, .keep_all = TRUE) %>%
  dplyr::select(any_of(c(
    "VariationID", "ClinicalSignificance",
    "LastEvaluated", "Description", "SubmittedPhenotypeInfo",
    "ReportedPhenotypeInfo", "ReviewStatus",
    "SubmittedGeneSymbol", "GeneSymbol", "vcf_id",
    "Origin", "OriginSimple", "Stars"
  )))

# extract all variants with conflicting interpretations
clinvar_conflicting <- variant_summary_df %>%
  dplyr::filter(ReviewStatus %in% c(
    "criteria provided, conflicting interpretations",
    "criteria provided, conflicting classifications"
  )) %>%
  pull(VariationID)

# extract all variants with >= P/LP call
plp_submissions <- submission_summary_df %>%
  dplyr::filter(ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic", "Pathogenic, low penetrance", "Likely pathogenic, low penetrance", "Established risk allele", "Likely risk allele")) %>%
  pull(VariationID)

# add column `clinvar_flag`, and make note of variants with conflicting interpretations and that have >= 1 P/LP call
submission_final_df <- submission_final_df %>%
  dplyr::mutate(clinvar_flag = case_when(
    VariationID %in% clinvar_conflicting & VariationID %in% plp_submissions ~ "Variant has conflicting interpretations in ClinVar and at least 1 submission with P/LP call",
    TRUE ~ "None"
  ))

table(submission_final_df$ClinicalSignificance)

write_tsv(
  submission_final_df,
  file.path(results_dir, glue::glue("resolved-clinvar-{date_submission}{concept_label}.tsv"))
)
