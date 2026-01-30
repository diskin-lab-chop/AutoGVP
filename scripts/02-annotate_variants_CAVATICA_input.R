################################################################################
# 02-annotate_variants_CAVATICA_input.R
# written by Ammar Naqvi & refactored by Saksham Phul
# updated 01/2026 by Patricia Sullivan
#
# This script annotates variants based on ClinVar and integrates a modified
# version of InterVar that involves adjustments of calls based on ACMG-AMP
# guidelines
#
# usage: Rscript 02-annotate_variants_CAVATICA_input.R --vcf <vcf file>
#                                       --clinvar  <ClinVar-selected-submissions.tsv>
#                                       --multianno <multianno file>
#                                       --intervar <intervar file>
#                                       --autopvs1 <autopvs1 file>
#                                       --output <string>
#                                       --outdir <output directory>
#                                       --sample_id <string>
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# parse parameters
option_list <- list(
  make_option(c("--vcf"),
    type = "character",
    help = "Input vcf file with VEP annotations"
  ),
  make_option(c("--intervar"),
    type = "character",
    help = "input intervar file"
  ),
  make_option(c("--multianno"),
    type = "character",
    help = "input multianno file"
  ),
  make_option(c("--autopvs1"),
    type = "character",
    help = "input autopvs1 file"
  ),
  make_option(c("--clinvar"),
    type = "character",
    help = "specific clinVar file (format: clinvar_20211225.vcf.gz)"
  ), ## https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022clinvar_20211225.vcf.gz
  make_option(c("--variant_summary"),
    type = "character",
    help = "variant_summary file (format: variant_summary_2023-02.txt)"
  ),
  make_option(c("--output"),
    type = "character", default = "out",
    help = "output name"
  ),
  make_option(c("--outdir"),
    type = "character", default = "../results",
    help = "output directory"
  ),
  make_option(c("--sample_id"),
    type = "character",
    help = "input sample bioassay id"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (read)
input_vcf_file <- opt$vcf
input_intervar_file <- opt$intervar
input_autopvs1_file <- opt$autopvs1
input_multianno_file <- opt$multianno
input_clinVar_file <- opt$clinvar
input_variant_summary <- opt$variant_summary
summary_level <- opt$summary_level_vcf
output_name <- opt$output
results_dir <- opt$outdir
sample_name <- opt$sample_id

# create results directory if it does not exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

## output files
output_tab_abr_file <- paste0(output_name, ".cavatica_input.annotations_report.abridged.tsv")

## allocate more memory capacity
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

input_vcf_file <- file.path("../data/test_VEP.vcf")
# Open vcf and read lines until a line without '#' is found
con <- file(input_vcf_file, "r")
skip_lines <- 0
while (grepl("^#", readLines(con, n = 1))) {
  skip_lines <- skip_lines + 1
}

## retrieve and store input vcf into table
vcf_df <- vroom(input_vcf_file, skip = skip_lines, delim = "\t", col_names = c("CHROM", "START", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"), show_col_types = FALSE) %>%
  dplyr::mutate(
    vcf_id = str_remove_all(paste(CHROM, "-", START, "-", REF, "-", ALT), " "),
    vcf_id = str_replace(vcf_id, "chr", "")
  )

input_variant_summary <- file.path("../results/ClinVar-selected-submissions.tsv")
## load in selected ClinVar submissions
clinvar_df <- vroom(input_variant_summary, show_col_types = FALSE) %>%
  dplyr::filter(vcf_id %in% vcf_df$vcf_id)

input_multianno_file <- file.path("../data/test_VEP.vcf.hg38_multianno.txt")
## get multianno file to add correct vcf_id in intervar table
multianno_df <- vroom(input_multianno_file, delim = "\t", trim_ws = TRUE, col_names = TRUE, show_col_types = FALSE) %>%
  dplyr::select(
    -End,
    -contains(c(
      "AF",
      "gnomad", "CLN",
      "score", "pred", "CADD", "Eigen",
      "100way", "30way", "GTEx"
    ))
  ) %>%
  mutate(
    vcf_id = str_remove_all(paste(Chr, "-", Otherinfo5, "-", Otherinfo7, "-", Otherinfo8), " "),
    vcf_id = str_replace(vcf_id, "chr", ""),
    var_id = str_remove_all(paste(Chr, "-", Start, "-", Ref, "-", Alt), " "),
    var_id = str_replace(var_id, "chr", "")
  ) %>%
  # remove coordiante, Otherinfo, gnomad, and clinVar-related columns
  dplyr::select(
    -Chr, -Ref, -Alt,
    -contains(("Otherinfo"))
  )

input_intervar_file <- file.path("../data/test_VEP.hg38_multianno.txt.intervar")
## add intervar table
intervar_df <- vroom(input_intervar_file, delim = "\t", trim_ws = TRUE, col_names = TRUE, show_col_types = FALSE) %>%
  dplyr::select(
    -`clinvar: Clinvar`,
    -contains(c("gnomad", "CADD", "Freq", "SCORE", "score", "ORPHA", "MIM", "rmsk", "GERP", "phylo"))
  ) %>%
  dplyr::mutate(var_id = paste0(`#Chr`, "-", Start, "-", Ref, "-", Alt)) %>%
  distinct(var_id, .keep_all = T) %>%
  # remove coordiante, Otherinfo, gnomad, and clinVar-related columns
  dplyr::select(
    -`#Chr`, -Start, -End, -Alt, -Ref
  ) %>%
  dplyr::select(any_of(c("Ref.Gene", "InterVar: InterVar and Evidence", "var_id"))) %>%
  ## add column for individual scores that will be re-calculated if we need to adjust using autoPVS1 result
  dplyr::mutate(
    evidencePVS1 = str_match(`InterVar: InterVar and Evidence`, "PVS1\\=(\\d+)\\s")[, 2],
    evidenceBA1 = str_match(`InterVar: InterVar and Evidence`, "BA1\\=(\\d+)\\s")[, 2],
    evidencePS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ","))))),
    evidencePM = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPM\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ","))))),
    evidencePP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))[-5])),
    evidenceBS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ","))))),
    evidenceBP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))[-6]))
    ## note: ignore PP5 score and BP6 score
  )

# combine intervar and multianno by var_id
intervar_multianno_df <- intervar_df %>%
  left_join(multianno_df, by = "var_id")

## combine the vcf, clinvar, intervar and multianno tables by the appropriate vcf id
clinvar_intervar_vcf_df <- vcf_df %>%
  left_join(intervar_multianno_df, by = "vcf_id") %>%
  left_join(clinvar_df %>% dplyr::select(any_of(c("vcf_id", "VariationID", "ClinicalSignificance", "ReviewStatus", "LastEvaluated", "clinvar_flag", "Origin", "OriginSimple"))), by = "vcf_id")


input_autopvs1_file <- file.path("../data/test_pbta.autopvs1.tsv")
## autopvs1 results
autopvs1_results <- vroom(input_autopvs1_file, col_names = TRUE, show_col_types = FALSE) %>%
  dplyr::select(vcf_id, Feature, criterion) %>%
  mutate(
    vcf_id = str_remove_all(vcf_id, " "),
    vcf_id = str_replace_all(vcf_id, "chr", "")
  ) %>%
  dplyr::filter(vcf_id %in% clinvar_intervar_vcf_df$vcf_id)


## merge autopvs1_results with vcf data, and filter for those variants that need intervar run
combined_tab_with_vcf_intervar <- autopvs1_results %>%
  inner_join(clinvar_intervar_vcf_df, by = "vcf_id") %>%
  ## indicate if recalculated
  dplyr::mutate(intervar_adjusted = if_else((evidencePVS1 == 0), "No", "Yes")) %>%
  dplyr::mutate(
    ## criteria to check intervar/autopvs1 to re-calculate and create a score column that will inform the new re-calculated final call
    # if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
    evidencePVS1 = if_else((criterion == "NF1" | criterion == "SS1" |
      criterion == "DEL1" | criterion == "DEL2" |
      criterion == "DUP1" | criterion == "IC1") & evidencePVS1 == 1, "1", evidencePVS1),

    # if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL4|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
    evidencePS = if_else((criterion == "NF3" | criterion == "NF5" |
      criterion == "SS3" | criterion == "SS5" |
      criterion == "SS8" | criterion == "SS10" | criterion == "DEL4" |
      criterion == "DEL8" | criterion == "DEL6" |
      criterion == "DEL10" | criterion == "DUP3" |
      criterion == "IC2") & evidencePVS1 == 1, as.numeric(evidencePS) + 1, as.double(evidencePS)),
    evidencePVS1 = if_else((criterion == "NF3" | criterion == "NF5" |
      criterion == "SS3" | criterion == "SS5" |
      criterion == "SS8" | criterion == "SS10" | criterion == "DEL4" |
      criterion == "DEL8" | criterion == "DEL6" |
      criterion == "DEL10" | criterion == "DUP3" |
      criterion == "IC2") & evidencePVS1 == 1, "0", evidencePVS1),

    # if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
    evidencePM = if_else((criterion == "NF6" | criterion == "SS6" |
      criterion == "SS9" | criterion == "DEL7" |
      criterion == "DEL11" | criterion == "IC3") & evidencePVS1 == 1, as.numeric(evidencePM) + 1, as.double(evidencePM)),
    evidencePVS1 = if_else((criterion == "NF6" | criterion == "SS6" |
      criterion == "SS9" | criterion == "DEL7" |
      criterion == "DEL11" | criterion == "IC3") & evidencePVS1 == 1, "0", evidencePVS1),

    # if criterion is IC4 then PVS1 = 0; PP = PP+1;
    evidencePP = if_else((criterion == "IC4") & evidencePVS1 == 1, as.numeric(evidencePP) + 1, as.double(evidencePP)),
    evidencePVS1 = if_else((criterion == "IC4") & evidencePVS1 == 1, "0", evidencePVS1),

    # if criterion is na|NF0|NF2|NF4|SS2|SS4|SS7|DEL3|DEL5|DEL9|DUP2|DUP4|DUP5|IC5 then PVS1 = 0;
    evidencePVS1 = if_else((criterion == "na" | criterion == "NF0" | criterion == "NF2" | criterion == "NF4" |
      criterion == "SS2" | criterion == "SS4" | criterion == "SS7" |
      criterion == "DEL3" | criterion == "DEL5" | criterion == "DEL9" |
      criterion == "DUP2" | criterion == "DUP4" | criterion == "DUP5" |
      criterion == "IC5") & evidencePVS1 == 1, 0, as.double(evidencePVS1)),

    ## adjust variables based on given rules described in README
    final_call_intervar = ifelse(intervar_adjusted == "No",
      sub(".*InterVar: ", "", sub("\\ P.*", "", `InterVar: InterVar and Evidence`)),
      ifelse((evidencePVS1 == 1) & (evidencePVS1 == 1 &
        ((evidencePS >= 1) |
          (evidencePM >= 2) |
          (evidencePM == 1 & evidencePP == 1) |
          (evidencePP >= 2)) &
        ((evidenceBA1) == 1 |
          (evidenceBS >= 2) |
          (evidenceBP >= 2) |
          (evidenceBS >= 1 & evidenceBP >= 1) |
          (evidenceBA1 == 1 & (evidenceBS >= 1 | evidenceBP >= 1)))), "Uncertain_significance",
      ifelse(((evidencePVS1 == 1) & (evidencePS >= 2) &
        ((evidenceBA1) == 1 |
          (evidenceBS >= 2) |
          (evidenceBP >= 2) |
          (evidenceBS >= 1 & evidenceBP >= 1) |
          (evidenceBA1 == 1 & (evidenceBS >= 1 | evidenceBP >= 1)))), "Uncertain_significance",
      ifelse(((evidencePVS1 == 1) & (evidencePS == 1 &
        (evidencePM >= 3 |
          (evidencePM == 2 & evidencePP >= 2) |
          (evidencePM == 1 & evidencePP >= 4))) &
        ((evidenceBA1) == 1 |
          (evidenceBS >= 2) |
          (evidenceBP >= 2) |
          (evidenceBS >= 1 & evidenceBP >= 1) |
          (evidenceBA1 == 1 & (evidenceBS >= 1 | evidenceBP >= 1)))), "Uncertain_significance",
      ifelse((((evidencePVS1 == 1) & (evidencePVS1 == 1 & evidencePM == 1) |
        (evidencePS == 1 & evidencePM >= 1) |
        (evidencePS == 1 & evidencePP >= 2) |
        (evidencePM >= 3) |
        (evidencePM == 2 & evidencePP >= 2) |
        (evidencePM == 1 & evidencePP >= 4)) &
        ((evidenceBA1) == 1 |
          (evidenceBS >= 2) |
          (evidenceBP >= 2) |
          (evidenceBS >= 1 & evidenceBP >= 1) |
          (evidenceBA1 == 1 & (evidenceBS >= 1 | evidenceBP >= 1)))), "Uncertain_significance",
      ifelse((evidencePVS1 == 1) & (evidencePVS1 == 1 &
        ((evidencePS >= 1) |
          (evidencePM >= 2) |
          (evidencePM == 1 & evidencePP == 1) |
          (evidencePP >= 2))), "Pathogenic",
      ifelse((evidencePVS1 == 1) & (evidencePS >= 2), "Pathogenic",
        ifelse((evidencePVS1 == 0) & (evidencePS == 1 &
          (evidencePM >= 3 |
            (evidencePM == 2 & evidencePP >= 2) |
            (evidencePM == 1 & evidencePP >= 4))), "Pathogenic",
        ifelse((evidencePVS1 == 1) & (evidencePVS1 == 1 & evidencePM == 1) |
          (evidencePS == 1 & evidencePM >= 1) |
          (evidencePS == 1 & evidencePP >= 2) |
          (evidencePM >= 3) |
          (evidencePM == 2 & evidencePP >= 2) |
          (evidencePM == 1 & evidencePP >= 4), "Likely_pathogenic",
        ifelse((evidencePVS1 == 1) & (evidenceBA1 == 1) |
          (evidenceBS >= 2), "Benign",
        ifelse((evidencePVS1 == 1) & (evidenceBS == 1 & evidenceBP == 1) |
          (evidenceBP >= 2), "Likely_benign", "Uncertain_significance")
        )
        )
        )
      )
      )
      )
      )
      )
      )
    )
  )


## merge tables together (clinvar and intervar) and write to file
master_tab <- clinvar_intervar_vcf_df %>%
  full_join(combined_tab_with_vcf_intervar[, grepl("vcf_id|intervar_adjusted|evidence|InterVar:|final_call_intervar", names(combined_tab_with_vcf_intervar))], by = "vcf_id") %>%
  left_join(autopvs1_results, by = "vcf_id") %>%
  # Make calls
  dplyr::mutate(
    intervar_adjusted = coalesce(intervar_adjusted, "No"),
    evidencePVS1 = coalesce(evidencePVS1.y, as.double(evidencePVS1.x)),
    evidenceBA1 = coalesce(as.double(evidenceBA1.y), as.double(evidenceBA1.x)),
    evidencePS = coalesce(as.double(evidencePS.y), as.double(evidencePS.x)),
    evidencePM = coalesce(as.double(evidencePM.y), as.double(evidencePM.x)),
    evidencePP = coalesce(as.double(evidencePP.y), as.double(evidencePP.x)),
    evidenceBS = coalesce(as.double(evidenceBS.y), as.double(evidenceBS.x)),
    evidenceBP = coalesce(as.double(evidenceBP.y), as.double(evidenceBP.x)),
    Intervar_evidence = coalesce(`InterVar: InterVar and Evidence.x`, `InterVar: InterVar and Evidence.y`),
  ) %>%
  dplyr::select(-c(
    evidencePVS1.x, evidencePVS1.y, evidenceBA1.x, evidenceBA1.y, evidencePS.x, evidencePS.y, evidencePM.x, evidencePM.y, evidencePP.x, evidencePP.y, evidenceBS.x, evidenceBS.y, evidenceBP.x, evidenceBP.y,
    `InterVar: InterVar and Evidence.x`, `InterVar: InterVar and Evidence.y`
  )) %>%
  dplyr::mutate(
    final_call_clinvar = ClinicalSignificance,
    final_call = coalesce(final_call_clinvar, final_call_intervar)
  )

## fix spelling and nomenclature inconsistencies
master_tab <- master_tab %>%
  dplyr::mutate(
    final_call = replace(final_call, final_call == "Likely benign", "Likely_benign"),
    final_call = replace(final_call, final_call == "Uncertain significance", "Uncertain_significance"),
    final_call = replace(final_call, final_call == "Benign PVS1", "Benign"),
    final_call = replace(final_call, final_call == "Pathogenic PVS1", "Pathogenic"),
    final_call = replace(final_call, final_call == "Likely pathogenic", "Likely_pathogenic")
  ) %>%
  distinct()


sample_name <- "Test"
## add column indicating final call source
master_tab <- master_tab %>%
  dplyr::mutate(Reasoning_for_call = case_when(
    is.na(final_call_clinvar) ~ "InterVar",
    TRUE ~ "ClinVar"
  )) %>%
  # modify `ClinVar_ClinicalSignificance` to equal `final_call` for ClinVar calls
  dplyr::mutate(ClinVar_ClinicalSignificance = case_when(
    Reasoning_for_call == "ClinVar" ~ final_call,
    TRUE ~ str_replace(final_call_clinvar, " ", "_")
  )) %>%
  dplyr::mutate(sample_id = sample_name) %>%
  dplyr::relocate(
    any_of(c(
      "sample_id", "CHROM", "START", "POS", "ID", "REF", "ALT",
      "final_call", "Reasoning_for_call",
      "Stars", "ClinVar_ClinicalSignificance", "Intervar_evidence"
    ))
  )

# write out to file
master_tab %>%
  write_tsv(
    file.path(results_dir, output_tab_abr_file),
    append = FALSE,
    quote = "none",
    col_names = TRUE
  )
