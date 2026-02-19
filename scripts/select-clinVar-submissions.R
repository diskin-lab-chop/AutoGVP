################################################################################
# select-clinVar-submissions.R
# written by Ryan Corbett, Patricia Sullivan
#
# This script selects unique ClinVar variant submission calls based on a list of
# predetermined criteria, to be used in AutoGVP for ClinVar variants that need
# resolving due to conflicting calls
#
# usage: select-clinVar-submissions.R --variant_summary <variant file>
#                                     --submission_summary <submission file>
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
    help = "variant_summary file (format: variant_summary.txt.gz)"
  ),
  make_option(c("--submission_summary"),
    type = "character",
    help = "specific submission summary file (format: submission_summary.txt.gz)"
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

input_submission_file <- file.path("../data/submission_summary_20260104.txt.gz")
input_variant_summary <- file.path("../data/variant_summary_20260104.txt.gz")
results_dir <- file.path("../refs/")
conceptID_file <- file.path("../refs/clinvar_cancer_concept_ids_20260130.txt")
conflict_res <- "latest"

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
      TRUE ~ NA
    ),
    Stars = case_when(
      ReviewStatus == "practice guideline" ~ 4,
      ReviewStatus == "reviewed by expert panel" ~ 3,
      ReviewStatus == "criteria provided, multiple submitters, no conflicts" ~ 2,
      ReviewStatus %in% c("criteria provided, conflicting classifications", "criteria provided, single submitter") ~ 1,
      TRUE ~ NA
    )
  )


# Load ClinVar submission summary file, which reports all submissions for each ClinVar variant

# Open the file and read lines until a line without '#' is found
con <- file(input_submission_file, "r")
skip_lines <- -1 # keep last line for column names
while (grepl("^#", readLines(con, n = 1))) {
  skip_lines <- skip_lines + 1
}

# Load submission file while skipping number of lines determined above
submission_summary_df <- vroom(input_submission_file,
  skip = skip_lines,
  delim = "\t",
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
  # Redefine `DateLastEvaluated`
  dplyr::mutate(
    DateLastEvaluated = case_when(
      DateLastEvaluated == "-" ~ NA_character_,
      TRUE ~ DateLastEvaluated
    ),
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
      "Benign",
      "Likely benign",
      "Uncertain significance"
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
  dplyr::mutate(LastEvaluated = coalesce(LastEvaluated_sub, LastEvaluated_var))

# Extract submissions that have no conflicts
variants_resolved <- submission_merged_df %>%
  filter(ClinSig_resolved != "Conflicting classifications of pathogenicity") %>%
  dplyr::arrange(desc(mdy(LastEvaluated))) %>%
  distinct(VariationID, .keep_all = T) %>%
  mutate(ClinSig_report = ClinSig_resolved) %>%
  arrange(VariationID)


# IF list of concept IDs provided -- filter remaining submissions to only those associated with concept IDs, and resolve conflicts by consensus, latest date, or severity
if (!is.null(conceptID_file)) {
  # read in concept IDs
  conceptIDs <- read_lines(conceptID_file)

  # Filter variants for those associated with concept IDs
  variants_with_conceptIDs <- submission_merged_df %>%
    dplyr::filter(ClinSig_resolved == "Conflicting classifications of pathogenicity") %>%
    dplyr::mutate(conceptID = unlist(lapply(strsplit(ReportedPhenotypeInfo, ":"), function(x) x[[1]]))) %>%
    dplyr::filter(conceptID %in% conceptIDs) %>%
    dplyr::select(-conceptID)
  
  length(unique(variants_with_conceptIDs$VariationID))

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

  # Identify variants with consensus calls
  consensus_calls_conceptIDs <- variants_with_conceptIDs %>%
    filter(!VariationID %in% variants_no_conflicts_conceptID) %>%
    # Here, we will group B+LB and P+LP together to identify majority calls
    dplyr::mutate(ClinSig = case_when(
      ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic") ~ "Pathogenic/Likely pathogenic",
      ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "Benign/Likely benign",
      ClinicalSignificance_sub == "Uncertain significance" ~ "Uncertain significance",
      TRUE ~ NA_character_
    )) %>%
    count(VariationID, ClinSig) %>%
    group_by(VariationID) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::filter(!VariationID %in% VariationID[duplicated(VariationID)])

  # extract variants resolved through consensus calling
  variants_resolved <- variants_with_conceptIDs %>%
    dplyr::mutate(ClinSig = case_when(
      ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic") ~ "Pathogenic/Likely pathogenic",
      ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "Benign/Likely benign",
      ClinicalSignificance_sub == "Uncertain significance" ~ "Uncertain significance",
      TRUE ~ NA_character_
    )) %>%
    dplyr::filter(glue::glue("{VariationID}-{ClinSig}") %in% glue::glue("{consensus_calls_conceptIDs$VariationID}-{consensus_calls_conceptIDs$ClinSig}")) %>%
    dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
    distinct(VariationID, .keep_all = T) %>%
    arrange(VariationID) %>%
    dplyr::mutate(ReviewStatus_var = glue::glue("{ReviewStatus_var}, submissions associated with conceptIDs, consensus call taken")) %>%
    dplyr::rename(ClinSig_report = ClinSig) %>%
    bind_rows(variants_resolved)

  # Resolve remaining conflicts for variants associated with concept IDs, either by taking date last evaluated or most severe call
  if (conflict_res == "latest") {
    variants_resolved <- variants_with_conceptIDs %>%
      dplyr::filter(!VariationID %in% consensus_calls_conceptIDs,
                    !VariationID %in% variants_no_conflicts_conceptID) %>%
      dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
      distinct(VariationID, .keep_all = T) %>%
      dplyr::mutate(
        ClinSig_report = ClinicalSignificance_sub,
        ReviewStatus_var = glue::glue("{ReviewStatus_var}, submissions associated with conceptIDs, latest date evaluated taken")) %>%
      bind_rows(variants_resolved)
  }

  if (conflict_res == "most_severe") {
    variants_resolved <- variants_with_conceptIDs %>%
      dplyr::filter(!VariationID %in% variants_consensus_call_conceptIDs$VariationID) %>%
      dplyr::mutate(ClinicalSignificance_sub = fct_relevel(
        ClinicalSignificance_sub,
        c(
          "Pathogenic", "Likely pathogenic",
          "Uncertain significance", 
          "Likely benign", "Benign"
        )
      )) %>%
      dplyr::arrange(VariationID, ClinicalSignificance_sub, desc(mdy(LastEvaluated_sub))) %>%
      distinct(VariationID, .keep_all = T) %>%
      dplyr::mutate(
        ClinSig_report = ClinicalSignificance_sub,
        ReviewStatus_var = glue::glue("{ReviewStatus_var}, submissions associated with conceptIDs, most severe call at latest date taken")) %>%
      bind_rows(variants_resolved)
  }
}

# Extract remaining variants with conflicts
conflicting_variants <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% variants_resolved$VariationID)

# Identify cases where a majority pathogenicity call has been made for variants
consensus_calls <- conflicting_variants %>%
  # Here, we will group B+LB and P+LP together to identify majority calls
  dplyr::mutate(ClinSig = case_when(
    ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic") ~ "Pathogenic/Likely pathogenic",
    ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "Benign/Likely benign",
    ClinicalSignificance_sub == "Uncertain significance" ~ "Uncertain significance",
    TRUE ~ NA_character_
  )) %>%
  count(VariationID, ClinSig) %>%
  group_by(VariationID) %>%
  dplyr::filter(n == max(n)) %>%
  dplyr::filter(!VariationID %in% VariationID[duplicated(VariationID)])

# Extract variants with consensus calls
variants_consensus_call <- conflicting_variants %>%
  dplyr::mutate(ClinSig = case_when(
    ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic") ~ "Pathogenic/Likely pathogenic",
    ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "Benign/Likely benign",
    ClinicalSignificance_sub == "Uncertain significance" ~ "Uncertain significance",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(glue::glue("{VariationID}-{ClinSig}") %in% glue::glue("{consensus_calls$VariationID}-{consensus_calls$ClinSig}")) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  distinct(VariationID, .keep_all = T) %>%
  arrange(VariationID) %>%
  dplyr::mutate(ReviewStatus_var = glue::glue("{ReviewStatus_var}, consensus call taken")) %>%
  dplyr::rename(ClinSig_report = ClinSig)

# Resolve conflicts in remaining variants by taking call at date last evaluated
variants_conflicts_latest <- conflicting_variants %>%
  dplyr::filter(!VariationID %in% variants_consensus_call) %>%
  dplyr::filter(ClinicalSignificance_var != ClinicalSignificance_sub) %>%
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
    ClinicalSignificance = ClinSig_report,
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
  distinct(vcf_id, .keep_all = T) %>%
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
  dplyr::filter(ClinicalSignificance %in% c("Pathogenic", "Likely pathogenic")) %>%
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
  file.path(results_dir, "ClinVar-selected-submissions.tsv")
)
