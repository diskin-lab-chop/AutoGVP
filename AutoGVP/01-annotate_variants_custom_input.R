################################################################################
# 01-annotate_variants_custom_input.R
# written by Ammar Naqvi & refactored by Saksham Phul
#
# This script annotates variants based on clinVar and integrates a modified
# version of InterVar that involves adjustments of calls based on ACMG-AMP
# guidelines
#
# usage: Rscript 01-annotate_variants_custom_input.R --vcf <VEP annotated vcf file>
#                                       --intervar <intervar file>
#                                       --autopvs1 <autopvs1 file>
#                                       --multianno <multianno file>
#                                       --clinvar  <e.g. clinvar_20211225.vcf.gz>
#                                       --variant_summary <variant_summary file>
#                                       --submission_summary <submission_summary file>
#                                       --output <string>
################################################################################

## load libraries
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
results_dir <- file.path(root_dir, "results")

# create results directory if it does not exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

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
    help = "specific clinVar file (format: clinvar_yyyymmdd.vcf.gz)"
  ), ## https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022clinvar_20211225.vcf.gz
  make_option(c("--variant_summary"),
    type = "character",
    help = "variant_summary  file (format: variant_summary_2023-02.txt)"
  ),
  make_option(c("--summary_level_vcf"),
    type = "character", default = "F",
    help = "summary_level T/F"
  ),
  make_option(c("--output"),
    type = "character", default = "out",
    help = "output name"
  )
)


opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (read)
input_clinVar_file <- opt$clinvar
input_intervar_file <- opt$intervar
input_autopvs1_file <- opt$autopvs1
input_vcf_file <- opt$vcf
input_variant_summary <- opt$variant_summary
input_multianno_file <- opt$multianno
summary_level <- opt$summary_level
output_name <- opt$output


## output files
output_tab_abr_file <- paste0(output_name, ".custom_input.annotations_report.abridged.tsv")

## allocate more memory capacity
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

address_ambiguous_calls <- function(results_tab_abridged) { ## address ambiguous calls (non L/LB/P/LP/VUS) by taking the InterVar final call

  results_tab_abridged <- results_tab_abridged %>%
    dplyr::mutate(new_call = case_when(
      is.na(final_call) | (final_call != "Pathogenic" &
        final_call != "Likely_benign" & final_call != "Likely_pathogenic" &
        final_call != "Uncertain_significance" & final_call != "Benign" &
        final_call != "Uncertain significance" & final_call != "Likely benign") ~ str_match(Intervar_evidence, "InterVar\\:\\s(\\w+\\s\\w+)*")[, 2],
      TRUE ~ NA_character_
    )) %>%
    dplyr::mutate(final_call = case_when(
      !is.na(new_call) ~ new_call,
      TRUE ~ final_call
    )) %>%
    dplyr::select(-new_call)

  return(results_tab_abridged)
}

address_conflicting_intrep <- function(clinvar_anno_vcf_df) { ## if conflicting intrep. take the call with most calls in CLNSIGCONF field
  for (i in 1:nrow(clinvar_anno_vcf_df)) {
    entry <- clinvar_anno_vcf_df[i, ]
    if (entry$Stars != "1NR") {
      next
    }

    conf_section <- str_match(entry$INFO, "CLNSIGCONF\\=.+\\;CLNVC") ## part to parse and count calls
    call_names <- c("Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign", "Uncertain_significance")

    P <- (str_match(conf_section, "Pathogenic\\((\\d+)\\)")[, 2])
    LP <- (str_match(conf_section, "Likely_pathogenic\\((\\d+)\\)")[, 2])
    B <- (str_match(conf_section, "Benign\\((\\d+)\\)")[, 2])
    LB <- (str_match(conf_section, "Likely_benign\\((\\d+)\\)")[, 2])
    U <- (str_match(conf_section, "Uncertain_significance\\((\\d+)\\)")[, 2])

    ## make vector out of possible calls to get max
    calls <- c(P, LP, B, LB, U)

    if (length(which(calls == max(calls, na.rm = TRUE))) > 1) {
      next
    }

    highest_ind <- which.max(calls)
    consensus_call <- call_names[highest_ind]

    clinvar_anno_vcf_df[i, ]$final_call <- consensus_call
  }
  return(clinvar_anno_vcf_df)
}

## make vcf dataframe and add vcf_if column
vcf_df <- vroom(input_vcf_file, comment = "#", delim = "\t", col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"), trim_ws = TRUE, show_col_types = FALSE) %>%
  mutate(
    vcf_id = str_remove_all(paste(CHROM, "-", POS, "-", REF, "-", ALT), " "),
    vcf_id = str_replace_all(vcf_id, "chr", "")
  )

## add clinvar table to this (INFO)
clinvar_anno_vcf_df <- vroom(input_clinVar_file, comment = "#", delim = "\t", col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), trim_ws = TRUE, show_col_types = FALSE) %>%
  # add vcf id column
  mutate(vcf_id = str_remove_all(paste(CHROM, "-", POS, "-", REF, "-", ALT), " ")) %>%
  semi_join(vcf_df, by = "vcf_id") %>%
  dplyr::mutate(
    vcf_id = str_replace_all(vcf_id, "chr", ""),
    # add star annotations to clinVar results table based on filters // ## default version
    Stars = ifelse(grepl("CLNREVSTAT\\=criteria_provided,_single_submitter", INFO), "1",
      ifelse(grepl("CLNREVSTAT\\=criteria_provided,_multiple_submitters", INFO), "2",
        ifelse(grepl("CLNREVSTAT\\=reviewed_by_expert_panel", INFO), "3",
          ifelse(grepl("CLNREVSTAT\\=practice_guideline", INFO), "4",
            ifelse(grepl("CLNREVSTAT\\=criteria_provided,_conflicting_interpretations", INFO), "1NR", "0")
          )
        )
      )
    ),
    ## extract the calls and put in own column
    final_call = str_match(INFO, "CLNSIG\\=(\\w+)([\\|\\/]\\w+)*\\;")[, 2]
  )

# clinvar_anno_vcf_df <- address_conflicting_intrep(clinvar_anno_vcf_df)

## store variants without clinvar info
clinvar_anti_join_vcf_df <- anti_join(vcf_df, clinvar_anno_vcf_df, by = "vcf_id") %>%
  dplyr::mutate(
    vcf_id = str_replace_all(vcf_id, "chr", ""),
    CHROM = str_replace_all(CHROM, "chr", "")
  ) %>%
  dplyr::rename(rs_id = ID)

## get latest calls from variant and submission summary files
variant_summary_df <- vroom(input_variant_summary)

## filter only those variants that need consensus call and find  call in submission table
entries_for_cc <- filter(clinvar_anno_vcf_df, Stars == "1NR", final_call != "Benign", final_call != "Pathogenic", final_call != "Likely_benign", final_call != "Likely_pathogenic", final_call != "Uncertain_significance")

entries_for_cc_in_submission <- inner_join(variant_summary_df, entries_for_cc, by = "vcf_id") %>%
  dplyr::mutate(final_call = ClinicalSignificance) %>%
  dplyr::select(vcf_id, ClinicalSignificance, final_call)

## one Star cases that are “criteria_provided,_single_submitter” that do NOT have the B, LB, P, LP, VUS call must also go to intervar
## modified: any cases that do NOT have the B, LB, P, LP, VUS call must also go to intervar
additional_intervar_cases <- filter(clinvar_anno_vcf_df, final_call != "Benign", final_call != "Pathogenic", final_call != "Likely_benign", final_call != "Likely_pathogenic", final_call != "Uncertain_significance") %>%
  anti_join(entries_for_cc_in_submission, by = "vcf_id") %>%
  anti_join(clinvar_anti_join_vcf_df, by = "vcf_id")

clinvar_anti_join_vcf_df <- clinvar_anti_join_vcf_df %>% mutate(QUAL = as.character(QUAL))

## filter only those variant entries that need an InterVar run (No Star) and add the additional intervar cases from above
entries_for_intervar <- filter(clinvar_anno_vcf_df, Stars == "0", na.rm = TRUE) %>%
  bind_rows((additional_intervar_cases)) %>%
  bind_rows(clinvar_anti_join_vcf_df) %>%
  distinct()

## get vcf ids that need intervar run
vcf_to_run_intervar <- entries_for_intervar$vcf_id

## get multianno file to add by correct vcf_id
multianno_df <- vroom(input_multianno_file, delim = "\t", trim_ws = TRUE, col_names = TRUE, show_col_types = FALSE) %>%
  mutate(
    vcf_id = str_remove_all(paste(Chr, "-", Otherinfo5, "-", Otherinfo7, "-", Otherinfo8), " "),
    vcf_id = str_replace_all(vcf_id, "chr", "")
  ) %>%
  group_by(vcf_id) %>%
  arrange(Chr, Start) %>%
  filter(row_number() == 1) %>%
  # remove coordiante, Otherinfo, gnomad, and clinVar-related columns
  dplyr::select(
    -Chr, -Start, -End, -Alt, -Ref,
    -contains(c(
      "Otherinfo", "gnomad", "CLN",
      "score", "pred", "CADD", "Eigen",
      "100way", "30way", "GTEx"
    ))
  ) %>%
  ungroup()

## add intervar table
clinvar_anno_intervar_vcf_df <- vroom(input_intervar_file, delim = "\t", trim_ws = TRUE, col_names = TRUE, show_col_types = TRUE) %>%
  # slice(-1) %>%
  dplyr::mutate(var_id = str_remove_all(paste(`#Chr`, "-", Start, "-", End, "-", Ref, "-", Alt), " ")) %>%
  group_by(var_id) %>%
  arrange(`#Chr`, Start) %>%
  filter(row_number() == 1) %>%
  # remove coordiante, Otherinfo, gnomad, and clinVar-related columns
  dplyr::select(
    -`#Chr`, -Start, -End, -Alt, -Ref, -`clinvar: Clinvar`,
    -contains(c("gnomad", "CADD", "Freq", "SCORE", "score", "ORPHA", "MIM", "rmsk"))
  ) %>%
  ungroup()

## exit if the total number of variants differ in these two tables to ensure we annotate with the correct vcf so we can match back to clinVar and other tables
if (tally(multianno_df) != tally(clinvar_anno_intervar_vcf_df)) {
  stop("intervar and multianno files of diff lengths")
}

## combine the intervar and multianno tables by the appropriate vcf id
clinvar_anno_intervar_vcf_df <-
  dplyr::mutate(multianno_df, clinvar_anno_intervar_vcf_df) %>%
  dplyr::filter(vcf_id %in% vcf_df$vcf_id)
# dplyr::select(any_of(c(
#   "vcf_id", "InterVar: InterVar and Evidence",
#   "Gene.refGene", "Ref.Gene", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene",
#   "CLNSIG", "CLNREVSTAT"
# )))

## populate consensus call variants with invervar info
entries_for_cc_in_submission_w_intervar <- inner_join(clinvar_anno_intervar_vcf_df, entries_for_cc_in_submission, by = "vcf_id") %>%
  # dplyr::select(any_of(c(
  #   "vcf_id", "InterVar: InterVar and Evidence",
  #   "Gene.refGene", "Ref.Gene", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene",
  #   "CLNSIG", "CLNREVSTAT"
  # ))) %>%
  dplyr::rename("Intervar_evidence" = `InterVar: InterVar and Evidence`)

## remove variants that we found in the submission file that were 1NR for intervar adjustment
clinvar_anno_intervar_vcf_df <- clinvar_anno_intervar_vcf_df %>%
  ## add column for individual scores that will be re-calculated if we need to adjust using autoPVS1 result

  ## note: ignore PP5 score and BP6 score
  dplyr::mutate(
    evidencePVS1 = str_match(`InterVar: InterVar and Evidence`, "PVS1\\=(\\d+)\\s")[, 2],
    evidenceBA1 = str_match(`InterVar: InterVar and Evidence`, "BA1\\=(\\d+)\\s")[, 2],
    evidencePS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ","))))),
    evidencePM = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPM\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ","))))),
    evidencePP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))[-5])),
    evidenceBS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ","))))),
    evidenceBP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))[-6]))
  ) %>%
  ## merge dataframe with clinvar_anno_vcf_df above
  left_join(vcf_df, by = "vcf_id") %>%
  left_join(clinvar_anno_vcf_df[, c("vcf_id", "Stars", "final_call")], by = "vcf_id")


## add back variants not found in clinVar db
# clinvar_anno_intervar_vcf_df <- bind_rows(clinvar_anno_intervar_vcf_df, clinvar_anti_join_vcf_df)

## autopvs1 results
autopvs1_results <- read_tsv(input_autopvs1_file, col_names = TRUE) %>%
  mutate(
    vcf_id = str_remove_all(paste(vcf_id), " "),
    vcf_id = str_replace_all(vcf_id, "chr", "")
  ) %>%
  dplyr::filter(vcf_id %in% clinvar_anno_intervar_vcf_df$vcf_id)

combined_tab_with_vcf_intervar <- autopvs1_results %>%
  inner_join(clinvar_anno_intervar_vcf_df, by = "vcf_id") %>%
  dplyr::filter(vcf_id %in% entries_for_intervar$vcf_id & !vcf_id %in% entries_for_cc_in_submission$vcf_id) %>%
  # dplyr::filter(vcf_id %in% entries_for_intervar$vcf_id)

  # combined_tab_for_intervar_cc_removed <- anti_join(combined_tab_with_vcf_intervar, entries_for_cc_in_submission, by = "vcf_id") %>%
  ## indicate if recalculated
  dplyr::mutate(intervar_adjusted = if_else((evidencePVS1 == 0), "No", "Yes")) %>%
  dplyr::mutate(
    ## criteria to check intervar/autopvs1 to re-calculate and create a score column that will inform the new re-calculated final call
    # if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
    evidencePVS1 = if_else((criterion == "NF1" | criterion == "SS1" |
      criterion == "DEL1" | criterion == "DEL2" |
      criterion == "DUP1" | criterion == "IC1") & evidencePVS1 == 1, "1", evidencePVS1),

    # if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
    evidencePS = if_else((criterion == "NF3" | criterion == "NF5" |
      criterion == "SS3" | criterion == "SS5" |
      criterion == "SS8" | criterion == "SS10" |
      criterion == "DEL8" | criterion == "DEL6" |
      criterion == "DEL10" | criterion == "DUP3" |
      criterion == "IC2") & evidencePVS1 == 1, as.numeric(evidencePS) + 1, as.double(evidencePS)),
    evidencePVS1 = if_else((criterion == "NF3" | criterion == "NF5" |
      criterion == "SS3" | criterion == "SS5" |
      criterion == "SS8" | criterion == "SS10" |
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

    # if criterion is na then PVS1 = 0;
    evidencePVS1 = if_else((criterion == "na") & evidencePVS1 == 1, 0, as.double(evidencePVS1)),

    ## adjust variables based on given rules described in README
    final_call = ifelse((evidencePVS1 == 1 &
      ((evidencePS >= 1) |
        (evidencePM >= 2) |
        (evidencePM == 1 & evidencePP == 1) |
        (evidencePP >= 2))), "Pathogenic",
    ifelse((evidencePS >= 2), "Pathogenic",
      ifelse((evidencePS == 1 &
        (evidencePM >= 3 |
          (evidencePM == 2 & evidencePP >= 2) |
          (evidencePM == 1 & evidencePP >= 4))), "Pathogenic",
      ifelse((evidencePVS1 == 1 & evidencePM == 1) |
        (evidencePS == 1 & evidencePM >= 1) |
        (evidencePS == 1 & evidencePP >= 2) |
        (evidencePM >= 3) |
        (evidencePM == 2 & evidencePP >= 2) |
        (evidencePM == 1 & evidencePP >= 4), "Likely_pathogenic",
      ifelse((evidenceBA1 == 1) |
        (evidenceBS >= 2), "Benign",
      ifelse((evidenceBS == 1 & evidenceBP == 1) |
        (evidenceBP >= 2), "Likely_benign", "Uncertain_significance")
      )
      )
      )
    )
    )
  )


## merge tables together (clinvar and intervar) and write to file
master_tab <- clinvar_anno_intervar_vcf_df %>%
  left_join(combined_tab_with_vcf_intervar[, grepl("vcf_id|intervar_adjusted|evidence|InterVar:|criterion|final_call", names(combined_tab_with_vcf_intervar))], by = "vcf_id") %>%
  #  left_join(combined_tab_for_intervar_cc_removed[, grepl("vcf_id|intervar_adjusted|evidence|InterVar:|criterion|final_call", names(combined_tab_for_intervar_cc_removed))], by = "vcf_id") %>%
  left_join(submission_final_df, by = "vcf_id")

master_tab <- master_tab %>%
  dplyr::mutate(
    intervar_adjusted = coalesce(intervar_adjusted, "No"),
    evidencePVS1 = coalesce(as.double(evidencePVS1.x, evidencePVS1.y)),
    evidenceBA1 = coalesce(as.double(evidenceBA1.x, evidenceBA1.y)),
    evidencePS = coalesce(as.double(evidencePS.x, evidencePS.y)),
    evidencePM = coalesce(as.double(evidencePM.x, evidencePM.y)),
    evidencePP = coalesce(as.double(evidencePP.x, evidencePP.y)),
    evidenceBS = coalesce(as.double(evidenceBS.x, evidenceBS.y)),
    evidenceBP = coalesce(as.double(evidenceBP.x, evidenceBP.y)),
    Intervar_evidence = coalesce(`InterVar: InterVar and Evidence.x`, `InterVar: InterVar and Evidence.y`),
    # replace second final call with the first one because we did not use interVar results
    final_call.x = if_else(evidencePVS1 == 0 & Stars == "0", final_call.y, final_call.x)
  )

## combine final calls into one choosing the appropriate final call
master_tab <- master_tab %>%
  dplyr::mutate(final_call = coalesce(final_call.x, final_call.y))

## remove older columns
master_tab <- master_tab %>% dplyr::select(-c(
  evidencePVS1.x, evidencePVS1.y, evidenceBA1.x, evidenceBA1.y, evidencePS.x, evidencePS.y, evidencePM.x, evidencePM.y, evidencePP.x, evidencePP.y, evidenceBS.x, evidenceBS.y, evidenceBP.x, evidenceBP.y,
  `InterVar: InterVar and Evidence.x`, `InterVar: InterVar and Evidence.y`, final_call.x, final_call.y
))

## reformat columns
master_tab <- full_join(master_tab, entries_for_cc_in_submission, by = "vcf_id") %>%
  dplyr::mutate(final_call = coalesce(final_call.y, final_call.x)) %>%
  full_join(entries_for_cc_in_submission_w_intervar[c("vcf_id", "Intervar_evidence")], by = "vcf_id") %>%
  dplyr::mutate(
    Intervar_evidence = coalesce(Intervar_evidence.y, Intervar_evidence.x),
    ClinVar_ClinicalSignificance = coalesce(ClinicalSignificance.x, ClinicalSignificance.y)
  ) %>%
  dplyr::select(
    -final_call.x, -final_call.y,
    -Intervar_evidence.x, -Intervar_evidence.y,
    -ClinicalSignificance.x, -ClinicalSignificance.y
  )


## address ambiguous calls (non L/LB/P/LP/VUS) by taking the InterVar final call
master_tab <- address_ambiguous_calls(master_tab)

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

## add column indicating final call source
master_tab <- master_tab %>%
  dplyr::mutate(Reasoning_for_call = case_when(
    vcf_id %in% vcf_to_run_intervar ~ "InterVar",
    TRUE ~ "ClinVar"
  )) %>%
  dplyr::relocate(any_of(c(
    "CHROM", "POS", "START", "ID", "REF", "ALT",
    "final_call", "Reasoning_for_call",
    "Stars", "ClinVar_ClinicalSignificance", "Intervar_evidence"
  )))

# write out to file
master_tab %>%
  write_tsv(
    file.path(results_dir, output_tab_abr_file),
    append = FALSE,
    quote = "none",
    col_names = TRUE
  )
