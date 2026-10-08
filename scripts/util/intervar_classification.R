################################################################################
# intervar_classification.R
#
# Shared InterVar re-classification used by 02-annotate_variants.R and
# fix_prior_results.R.
#
# classify_intervar_call() takes summed ACMG evidence counts (PP5 and BP6
# already excluded from `PP` and `BP`) and returns the InterVar call.
# Pathogenic and benign evidence are evaluated separately; if both are met the
# evidence is conflicting and the call is Uncertain_significance (as in
# InterVar's `classify()`).
################################################################################

classify_intervar_call <- function(PVS1, PS, PM, PP, BA1, BS, BP) {
  path_p <- (PVS1 == 1 & (PS >= 1 | PM >= 2 | (PM == 1 & PP == 1) | PP >= 2)) |
    PS >= 2 |
    (PS == 1 & (PM >= 3 | (PM == 2 & PP >= 2) | (PM == 1 & PP >= 4)))

  path_lp <- (PVS1 == 1 & PM == 1) |
    (PS == 1 & PM >= 1) |
    (PS == 1 & PP >= 2) |
    PM >= 3 |
    (PM == 2 & PP >= 2) |
    (PM == 1 & PP >= 4)

  benign_b <- BA1 == 1 | BS >= 2
  benign_lb <- (BS >= 1 & BP >= 1) | BP >= 2

  dplyr::case_when(
    (path_p | path_lp) & (benign_b | benign_lb) ~ "Uncertain_significance",
    path_p ~ "Pathogenic",
    path_lp ~ "Likely_pathogenic",
    benign_b ~ "Benign",
    benign_lb ~ "Likely_benign",
    TRUE ~ "Uncertain_significance"
  )
}
