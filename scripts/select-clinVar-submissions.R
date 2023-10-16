################################################################################
# select-clinVar-submissions.R
# written by Ryan Corbett
#
# This script selects unique clinVar variant submission calls based on a list of
# predetermined criteria, to be used in AutoGVP for ClinVar variants that need
# resolving due to conflicting calls
#
# usage: select-clinVar-submissions.R --variant_summary <variant file>
#                                       --submission_summary <submission file>
#
# NOTE: this script must be run BEFORE running run_autogvp.sh
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path(root_dir, "data")


# parse parameters
option_list <- list(
  make_option(c("--variant_summary"),
    type = "character",
    help = "variant_summary file (format: variant_summary.txt.gz)"
  ),
  make_option(c("--submission_summary"),
    type = "character",
    help = "specific submission summary file (format: submission_summary.txt.gz)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (reqd)
input_submission_file <- opt$submission_summary
input_variant_summary <- opt$variant_summary


## load variant summary file, which reports latest clinVar consensus calls for each variant
variant_summary_df <- vroom(input_variant_summary,
  delim = "\t",
  col_types = c(ReferenceAlleleVCF = "c", AlternateAlleleVCF = "c", PositionVCF = "i", VariationID = "n"),
  show_col_types = FALSE
) %>%
  # retain only variants mapped to hg38
  dplyr::filter(Assembly == "GRCh38" & ReferenceAlleleVCF != "na" & AlternateAlleleVCF != "na") %>%
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
  dplyr::filter(!ReviewStatus %in% c("no assertion provided", "no assertion criteria provided"))

# Load clinVar submission summary file, which reports all submissions for each clinVar variant
submission_summary_df <- vroom(input_submission_file,
  comment = "#", delim = "\t",
  col_names = c(
    "VariationID", "ClinicalSignificance", "DateLastEvaluated",
    "Description", "SubmittedPhenotypeInfo", "ReportedPhenotypeInfo",
    "ReviewStatus", "CollectionMethod", "OriginCounts", "Submitter",
    "SCV", "SubmittedGeneSymbol", "ExplanationOfInterpretation"
  ),
  show_col_types = F
) %>%
  # Redefine `DateLastEvaluated`
  dplyr::mutate(DateLastEvaluated = case_when(
    DateLastEvaluated == "-" ~ NA_character_,
    TRUE ~ DateLastEvaluated
  )) %>%
  dplyr::filter(!ReviewStatus %in% c("no assertion provided", "no assertion criteria provided"))

# merge submission_summary and variant_summary info
submission_merged_df <- submission_summary_df %>%
  dplyr::rename("LastEvaluated" = DateLastEvaluated) %>%
  left_join(variant_summary_df,
    by = "VariationID",
    relationship = "many-to-many", suffix = c("_sub", "_var")
  ) %>%
  dplyr::mutate(LastEvaluated = coalesce(LastEvaluated_sub, LastEvaluated_var)) %>%
  dplyr::filter(!is.na(vcf_id))

# Extract submissions that match variant consensus call and are 2+ stars
variants_no_conflict_expert <- submission_merged_df %>%
  filter(ReviewStatus_var %in% c(
    "practice guideline", "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts"
  )) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()
#  filter(ClinicalSignificance_sub == ClinicalSignificance_var | is.na(ClinicalSignificance_var) | (grepl("Pathogenic|Likely pathogenic", ClinicalSignificance_sub) & grepl("Pathogenic|Likely pathogenic", ClinicalSignificance_var)) | (grepl("Benign|Likely benign", ClinicalSignificance_sub) & grepl("Benign|Likely benign", ClinicalSignificance_var)) | (grepl("Uncertain significance", ClinicalSignificance_sub) & grepl("Uncertain significance", ClinicalSignificance_var)))

# Identify VariationIDs with no ClinSig conflicts between variant and submission summary at date last evaluated
variants_no_conflicts <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% variants_no_conflict_expert$VariationID) %>%
  dplyr::filter(ClinicalSignificance_sub == ClinicalSignificance_var | is.na(ClinicalSignificance_var) | (grepl("Pathogenic|Likely pathogenic", ClinicalSignificance_sub) & grepl("Pathogenic|Likely pathogenic", ClinicalSignificance_var)) | (grepl("Benign|Likely benign", ClinicalSignificance_sub) & grepl("Benign|Likely benign", ClinicalSignificance_var)) | (grepl("Uncertain significance", ClinicalSignificance_sub) & grepl("Uncertain significance", ClinicalSignificance_var))) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()

# Identify cases where a majority pathogenicity call has been made for variants
consensus_calls <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% c(variants_no_conflict_expert$VariationID, variants_no_conflicts$VariationID)) %>%
  # Here, we will group B+LB and P+LP together to identify majority calls
  dplyr::mutate(ClinSig = case_when(
    ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic") ~ "P/LP",
    ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "B/LB",
    ClinicalSignificance_sub == "Uncertain significance" ~ "VUS",
    TRUE ~ NA_character_
  )) %>%
  count(VariationID, ClinSig) %>%
  group_by(VariationID) %>%
  dplyr::filter(n == max(n)) %>%
  dplyr::filter(!VariationID %in% VariationID[duplicated(VariationID)])

# Extract variants with majority calls
variants_consensus_call <- submission_merged_df %>%
  dplyr::mutate(ClinSig = case_when(
    ClinicalSignificance_sub %in% c("Pathogenic", "Likely pathogenic") ~ "P/LP",
    ClinicalSignificance_sub %in% c("Benign", "Likely benign") ~ "B/LB",
    ClinicalSignificance_sub == "Uncertain significance" ~ "VUS",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(glue::glue("{VariationID}-{ClinSig}") %in% glue::glue("{consensus_calls$VariationID}-{consensus_calls$ClinSig}")) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::mutate(ReviewStatus_var = glue::glue("{ReviewStatus_var}, consensus call taken")) %>%
  dplyr::select(-ClinSig)

# Identify variants with conflicting ClinSigs, and retain call at most recent date evaluated
variants_conflicts_latest <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% c(variants_no_conflict_expert$VariationID, variants_no_conflicts$VariationID, variants_consensus_call$VariationID)) %>%
  dplyr::filter(ClinicalSignificance_var != ClinicalSignificance_sub) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated_sub))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::mutate(ReviewStatus_var = glue::glue("{ReviewStatus_var}, latest date evaluated call taken"))

# create final df and take ClinSig calls from submission summary
submission_final_df <- variants_no_conflicts %>%
  bind_rows(variants_no_conflict_expert, variants_consensus_call, variants_conflicts_latest) %>%
  dplyr::mutate(
    ClinicalSignificance = ClinicalSignificance_sub,
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
    "SubmittedGeneSymbol", "GeneSymbol", "vcf_id"
  )))

write_tsv(
  submission_final_df,
  file.path(input_dir, "ClinVar-selected-submissions.tsv")
)
