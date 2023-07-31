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
# NOTE: this script must be run BEFORE running run-autogvp.sh
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
analysis_dir <- file.path(root_dir, "AutoGVP")
input_dir <- file.path(analysis_dir, "input")


# parse parameters
option_list <- list(
  make_option(c("--variant_summary"),
    type = "character",
    help = "variant_summary file (format: variant_summary_2023-02.txt)"
  ),
  make_option(c("--submission_summary"),
    type = "character",
    help = "specific submission summary file (format: submission_summary.txt)"
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
  )

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
)

# remove submissions with missing columns (will have NA SubmittedGeneSymbol)
submission_summary_df <- submission_summary_df[!is.na(submission_summary_df$SubmittedGeneSymbol), ]

# Redefine `DateLastEvaluated`
submission_summary_df <- submission_summary_df %>%
  dplyr::mutate(DateLastEvaluated = case_when(
    DateLastEvaluated == "-" ~ NA_character_,
    TRUE ~ DateLastEvaluated
  ))

# merge submission_summary and variant_summary info
submission_merged_df <- submission_summary_df %>%
  dplyr::rename("LastEvaluated" = DateLastEvaluated) %>%
  left_join(variant_summary_df, by = c("VariationID")) %>%
  dplyr::mutate(LastEvaluated = coalesce(LastEvaluated.x, LastEvaluated.y)) %>%
  dplyr::filter(!is.na(vcf_id))

# Extract submissions that match variant consensus call and are reviewed by expert panel
variants_no_conflict_expert <- submission_merged_df %>%
  filter(ReviewStatus.x == "reviewed by expert panel") %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated.x))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup() %>%
  filter(ClinicalSignificance.x == ClinicalSignificance.y | is.na(ClinicalSignificance.y) | (grepl("Pathogenic|Likely pathogenic", ClinicalSignificance.x) & grepl("Pathogenic|Likely pathogenic", ClinicalSignificance.y)) | (grepl("Benign|Likely benign", ClinicalSignificance.x) & grepl("Benign|Likely benign", ClinicalSignificance.y)) | (grepl("Uncertain significance", ClinicalSignificance.x) & grepl("Uncertain significance", ClinicalSignificance.y)))

# Identify VariationIDs with no ClinSig conflicts between variant and submission summary at date last evaluated
variants_no_conflicts <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% variants_no_conflict_expert$VariationID) %>%
  dplyr::filter(ClinicalSignificance.x == ClinicalSignificance.y | is.na(ClinicalSignificance.y) | (grepl("Pathogenic|Likely pathogenic", ClinicalSignificance.x) & grepl("Pathogenic|Likely pathogenic", ClinicalSignificance.y)) | (grepl("Benign|Likely benign", ClinicalSignificance.x) & grepl("Benign|Likely benign", ClinicalSignificance.y)) | (grepl("Uncertain significance", ClinicalSignificance.x) & grepl("Uncertain significance", ClinicalSignificance.y))) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated.x))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()

# Identify cases where a majority pathnogenicity calls has been made for variants
consensus_calls <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% c(variants_no_conflict_expert$VariationID, variants_no_conflicts$VariationID)) %>%
  count(VariationID, ClinicalSignificance.x) %>%
  group_by(VariationID) %>%
  dplyr::filter(n == max(n)) %>%
  dplyr::filter(!VariationID %in% VariationID[duplicated(VariationID)])

# Extract variants with majority calls
variants_consensus_call <- submission_merged_df %>%
  dplyr::filter(glue::glue("{VariationID}-{ClinicalSignificance.x}") %in% glue::glue("{consensus_calls$VariationID}-{consensus_calls$ClinicalSignificance.x}")) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated.x))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()

# Identify variants with conflicting ClinSigs, but where a P-LP call has an associated phenotypeInfo
variants_conflicts_phenoInfo <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% c(variants_no_conflict_expert$VariationID, variants_no_conflicts$VariationID, variants_consensus_call$VariationID)) %>%
  dplyr::filter(ClinicalSignificance.y != ClinicalSignificance.x) %>%
  dplyr::filter(grepl("Pathogenic|Likely pathogenic", ClinicalSignificance.x) & (!is.na(SubmittedPhenotypeInfo) | !is.na(ReportedPhenotypeInfo)) & !grepl("not provided|not specified", ReportedPhenotypeInfo)) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated.x))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()

# Identify variants with conflicting ClinSigs, and retain call at last date evaluated
variants_conflicts_latest <- submission_merged_df %>%
  dplyr::filter(!VariationID %in% c(variants_no_conflict_expert$VariationID, variants_no_conflicts$VariationID, variants_consensus_call$VariationID, variants_conflicts_phenoInfo$VariationID)) %>%
  dplyr::filter(ClinicalSignificance.y != ClinicalSignificance.x) %>%
  group_by(VariationID) %>%
  dplyr::arrange(desc(mdy(LastEvaluated.x))) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()

# create final df and take ClinSig calls from submission summary
submission_final_df <- variants_no_conflicts %>%
  bind_rows(variants_no_conflict_expert, variants_consensus_call, variants_conflicts_phenoInfo, variants_conflicts_latest) %>%
  dplyr::mutate(
    ClinicalSignificance = ClinicalSignificance.x,
    ReviewStatus = ReviewStatus.y,
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
    "VariationID", "ClinicalSignificance", "ClinicalSignificance",
    "LastEvaluated", "Description", "SubmittedPhenotypeInfo",
    "ReportedPhenotypeInfo", "ReviewStatus",
    "SubmittedGeneSymbol", "GeneSymbol", "vcf_id"
  )))

write_tsv(
  submission_final_df,
  file.path(input_dir, "ClinVar-selected-submissions.tsv")
)
