################################################################################
# update_intervar.R
# written by Ryan Corbett
#
# This script loads the specified intervar df and updates PS1 and PM5 criteria
# and intervar calls based on clinvar entries in supplied resolved clinvar
# interpretations
#
# usage: update_intervar.R --intervar_file <intervar file>
#                          --clinvar_file <resolved clinvar file>
#                          --clinvar_hgvs4_file <clinvar hgvs4 variation file>
#
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
  make_option(c("--intervar_file"),
    type = "character",
    help = "Input intervar file"
  ),
  make_option(c("--clinvar_file"),
    type = "character",
    help = "resolved clinvar interpretations"
  ),
  make_option(c("--clinvar_hgvs4_file"),
    type = "character",
    help = "clinvar hgvs4 variation file"
  ),
  make_option(c("--outdir"),
    type = "character", default = "../results",
    help = "output directory"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (read)
intervar_file <- opt$intervar_file
clinvar_file <- opt$clinvar_file
clinvar_hgvs4_file <- opt$clinvar_hgvs4_file
results_dir <- opt$outdir
output_file <- glue::glue("{basename(intervar_file)}_updated")

# Load InterVar output and ensure chromosome values are stored as character
# strings to avoid type mismatches during downstream joins
intervar_df <- read_tsv(intervar_file, show_col_types = FALSE) %>%
  dplyr::mutate(`#Chr` = as.character(`#Chr`))

# Parse genomic coordinates from the ClinVar vcf_id field
# (chr-pos-ref-alt format) into separate columns for variant matching
clinvar_df <- vroom::vroom(
  clinvar_file,
  show_col_types = FALSE
) %>%
  tidyr::separate_wider_delim(
    vcf_id,
    delim = "-",
    names = c(
      "chr_clinvar",
      "Start_clinvar",
      "Ref_clinvar",
      "Alt_clinvar"
    )
  ) %>%
  mutate(
    Start_clinvar = as.numeric(Start_clinvar)
  )

# function to convert clinvar HGVSp annotations to one-letter AA abbreviations
# to match intervar abbreviations
convert_hgvsp_to_one_letter <- function(hgvsp) {
  aa_map <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D",
    Cys = "C", Gln = "Q", Glu = "E", Gly = "G",
    His = "H", Ile = "I", Leu = "L", Lys = "K",
    Met = "M", Phe = "F", Pro = "P", Ser = "S",
    Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Ter = "*"
  )

  for (aa3 in names(aa_map)) {
    hgvsp <- str_replace_all(hgvsp, aa3, aa_map[[aa3]])
  }

  hgvsp
}

# subset intervar df for missense variants
# extract relevant criteria to separate columns
# append clinvar annotations
intervar_missense_df <- intervar_df %>%
  dplyr::filter(ExonicFunc.refGene == "nonsynonymous SNV") %>%
  # parse out AAChanges for each transcript in `AAChange.knownGene` into unique rows
  tidyr::separate_longer_delim(
    AAChange.knownGene,
    delim = ","
  ) %>%
  # Extract previously assigned InterVar PS1 and PM5 evidence values
  # so they can be compared against updated ClinVar-derived values
  dplyr::mutate(
    PS1_old = as.integer(str_match(`InterVar: InterVar and Evidence`, "PS=\\[([01])")[, 2]),
    PM5_old = as.integer(
      str_match(
        `InterVar: InterVar and Evidence`,
        "PM=\\[[01],\\s*[01],\\s*[01],\\s*[01],\\s*([01])"
      )[, 2]
    ),
    # Extract HGVS protein change annotations from AAChange.knownGene.
    # Following transcript expansion, each row should contain only a
    # single transcript annotation in the format:
    # Gene:Transcript:Exon:cDNA_change:Protein_change
    HGVSp = if_else(
      AAChange.knownGene == ".",
      NA_character_,
      vapply(
        strsplit(AAChange.knownGene, ":", fixed = TRUE),
        function(x) if (length(x) >= 5) x[5] else NA_character_,
        character(1)
      )
    )
  ) %>%
  dplyr::mutate(HGVSp = str_remove(HGVSp, ",")) %>%
  dplyr::select(
    `#Chr`, Start, End, Ref, Alt,
    HGVSp, `clinvar: Clinvar`,
    `InterVar: InterVar and Evidence`,
    PS1_old, PM5_old
  ) %>%
  # Match each InterVar variant to ClinVar records occurring at the
  # same genomic position. Multiple ClinVar records may be associated
  # with a single InterVar variant.
  left_join(
    clinvar_df %>%
      dplyr::select(
        chr_clinvar, Start_clinvar,
        Ref_clinvar, Alt_clinvar,
        VariationID,
        ClinicalSignificance
      ),
    by = c(
      "#Chr" = "chr_clinvar",
      "Start" = "Start_clinvar"
    ),
    relationship = "many-to-many"
  ) %>%
  # retain only SNVs (PS1 and PM5 only applied to missense variants)
  dplyr::filter(
    nchar(Ref_clinvar) == 1,
    nchar(Alt_clinvar) == 1
  ) %>%
  dplyr::rename(
    ClinVar_VariationID = VariationID,
    ClinicalSignificance_old_clinvar = `clinvar: Clinvar`,
    HGVSp = HGVSp,
    ClinicalSignificance_new_clinvar = ClinicalSignificance
  )

# Load ClinVar HGVS4Variation file, which provides protein-level HGVS
# annotations (ProteinChange) linked to ClinVar VariationIDs.
# This file is typically tens of GB uncompressed, but only rows whose
# VariationID matches an InterVar missense variant are actually needed.
# Pre-filter with awk (a fast, single C-level pass over the decompressed
# stream) so that only the small matching subset ever reaches vroom,
# instead of loading and indexing the entire file in R.
intervar_ids <- unique(intervar_missense_df$ClinVar_VariationID)

ids_file <- tempfile()
writeLines(as.character(intervar_ids), ids_file)

decompress_cmd <- if (nzchar(Sys.which("zcat"))) "zcat" else "gzcat"

hgvs4_cols <- c(
  "Symbol", "GeneID", "VariationID", "AlleleID", "Type", "Assembly",
  "NucleotideExpression", "NucleotideChange", "ProteinExpression",
  "ProteinChange", "UsedForNaming", "Submitted", "OnRefSeqGene"
)

# VariationID is column 3 of the HGVS4Variation file; keep only data
# lines (drop the leading "#..." comment/header lines) whose VariationID
# is in the InterVar-matched set
filter_cmd <- sprintf(
  "%s %s | awk -F'\\t' 'NR==FNR { gsub(/\\r$/, \"\"); ids[$1]; next } { gsub(/\\r$/, \"\") } !/^#/ && ($3 in ids)' %s -",
  decompress_cmd,
  shQuote(clinvar_hgvs4_file),
  shQuote(ids_file)
)

hgvs4_variation_df <- vroom::vroom(
  pipe(filter_cmd),
  delim = "\t",
  col_names = hgvs4_cols,
  col_select = c(VariationID, Assembly, ProteinChange),
  show_col_types = FALSE
)

unlink(ids_file)

# Verify that HGVS annotations were found for all ClinVar variants.
# Missing VariationIDs usually indicate an outdated HGVS4Variation file.
if (length(unique(intervar_missense_df$ClinVar_VariationID)) > length(unique(hgvs4_variation_df$VariationID))) {
  print("Warning: HGVS annotations not found for some ClinVar variants. Please ensure you are supplying a recent hgvs4variation.txt.gz file from ClinVar")
}

# Retain only protein-level annotations from the HGVS4Variation file.
# Assembly == "na" corresponds to protein annotations rather than
# genomic or transcript-level representations.
hgvs4_variation_df <- hgvs4_variation_df %>%
  # only retain entries with unique AA changes
  dplyr::filter(
    Assembly == "na",
    ProteinChange != "-"
  ) %>%
  distinct(VariationID, ProteinChange)

# Append ClinVar protein change annotations to each matched variant
# and convert amino acid abbreviations to match InterVar formatting
# (e.g. Trp507Arg -> W507R)
intervar_missense_df <- intervar_missense_df %>%
  left_join(
    hgvs4_variation_df %>%
      dplyr::select(VariationID, ProteinChange),
    by = c("ClinVar_VariationID" = "VariationID"),
    relationship = "many-to-many"
  ) %>%
  dplyr::rename(HGVSp_clinvar = ProteinChange) %>%
  # convert to one-letter AA abbreviations to match intervar
  dplyr::mutate(HGVSp_clinvar = convert_hgvsp_to_one_letter(HGVSp_clinvar))

# Helper function to extract amino acid residue positions from HGVS
# protein annotations (e.g. p.W343R -> 343)
extract_aa_pos <- function(x) {
  as.numeric(str_extract(x, "(?<=\\D)\\d+"))
}

# Calculate new PS1, PM5 criteria based on ClinVar entries
intervar_missense_df <- intervar_missense_df %>%
  mutate(
    # Get AA positions
    aa_pos = extract_aa_pos(HGVSp),
    aa_pos_clinvar = extract_aa_pos(HGVSp_clinvar),

    # determine if clinvar variant is PLP (PS1, PM5 criteria)
    clinvar_plp = str_detect(
      ClinicalSignificance_new_clinvar,
      regex("pathogenic|likely pathogenic", ignore_case = TRUE)
    ) &
      !str_detect(
        ClinicalSignificance_new_clinvar,
        regex("benign", ignore_case = TRUE)
      ),

    # determine if clinvar variant is a different nt change
    different_nt =
      Ref != Ref_clinvar |
        Alt != Alt_clinvar,

    # Compute relationships between the InterVar variant and matched
    # ClinVar variants that are required for ACMG PS1 and PM5 evaluation
    same_protein = HGVSp == HGVSp_clinvar,
    same_codon =
      aa_pos == aa_pos_clinvar
  ) %>%
  group_by(`#Chr`, Start, Ref, Alt, HGVSp) %>%
  mutate(
    # Evaluate ACMG evidence criteria for each InterVar-ClinVar match:
    # PS1: same amino acid substitution produced by a different nucleotide change
    # PM5: different amino acid substitution affecting the same codon
    PS1_new = as.integer(
      clinvar_plp &
        different_nt &
        same_protein
    ),

    # PM5: P/LP variant at same position, diff nt, diff AA change
    PM5_new = as.integer(
      clinvar_plp &
        different_nt &
        same_codon &
        !same_protein
    )
  ) %>%
  ungroup()

# Aggregate evidence across all ClinVar matches associated with the
# same InterVar variant. A single supporting ClinVar record is sufficient
# to activate PS1 or PM5 for the variant.
variant_summary <- intervar_missense_df %>%
  dplyr::mutate(
    PS1_support = clinvar_plp &
      different_nt &
      same_protein,
    PM5_support = clinvar_plp &
      different_nt &
      same_codon &
      !same_protein
  ) %>%
  group_by(`#Chr`, Start, Ref, Alt, HGVSp) %>%
  summarise(
    PS1_new = as.integer(any(PS1_support, na.rm = TRUE)),
    PM5_new = as.integer(any(PM5_support, na.rm = TRUE)),
    PS1_ClinVarIDs = paste(
      unique(ClinVar_VariationID[PS1_support]),
      collapse = ";"
    ),
    PM5_ClinVarIDs = paste(
      unique(ClinVar_VariationID[PM5_support]),
      collapse = ";"
    ),
    .groups = "drop"
  )

# Reduce back to one row per InterVar variant and merge the aggregated
# PS1/PM5 evidence assignments calculated from ClinVar
intervar_unique <- intervar_missense_df %>%
  dplyr::select(
    -PS1_new, -PM5_new,
    -Ref_clinvar, -Alt_clinvar,
    -ClinVar_VariationID, -ClinicalSignificance_new_clinvar,
    -HGVSp_clinvar, -aa_pos,
    -aa_pos_clinvar, -clinvar_plp,
    -different_nt, -same_protein, -same_codon
  ) %>%
  distinct(`#Chr`, Start, Ref, Alt, HGVSp, .keep_all = TRUE) %>%
  left_join(
    variant_summary,
    by = c("#Chr", "Start", "Ref", "Alt", "HGVSp")
  )

# function to update intervar to match formatting of `InterVar: InterVar and Evidence`
update_intervar <- function(intervar_string, PS1_new, PM5_new) {
  # Parse individual ACMG evidence codes from the InterVar annotation
  # string so that PS1 and PM5 can be updated and the classification
  # recalculated
  PVS1 <- as.numeric(str_match(intervar_string, "PVS1=(\\d)")[, 2])

  PS <- str_match(intervar_string, "PS=\\[([^]]+)\\]")[, 2] |>
    str_split(",\\s*") |>
    unlist() |>
    as.numeric()

  PM <- str_match(intervar_string, "PM=\\[([^]]+)\\]")[, 2] |>
    str_split(",\\s*") |>
    unlist() |>
    as.numeric()

  PP <- str_match(intervar_string, "PP=\\[([^]]+)\\]")[, 2] |>
    str_split(",\\s*") |>
    unlist() |>
    as.numeric()

  BA1 <- as.numeric(str_match(intervar_string, "BA1=(\\d)")[, 2])

  BS <- str_match(intervar_string, "BS=\\[([^]]+)\\]")[, 2] |>
    str_split(",\\s*") |>
    unlist() |>
    as.numeric()

  BP <- str_match(intervar_string, "BP=\\[([^]]+)\\]")[, 2] |>
    str_split(",\\s*") |>
    unlist() |>
    as.numeric()

  # save original values
  old_PS1 <- PS[1]
  old_PM5 <- PM[5]

  # update evidence
  PS[1] <- max(PS[1], PS1_new)
  PM[5] <- max(PM[5], PM5_new)

  # if nothing changed, return original string
  if (PS[1] == old_PS1 && PM[5] == old_PM5) {
    return(intervar_string)
  }

  nPS <- sum(PS)
  nPM <- sum(PM)
  nPP <- sum(PP)
  nBS <- sum(BS)
  nBP <- sum(BP)

  # Recalculate the InterVar classification using the updated ACMG
  # evidence profile and standard InterVar decision rules
  new_class <- case_when(
    BA1 == 1 ~ "Benign",
    nBS >= 2 ~ "Benign",
    (nBS == 1 & nBP >= 1) | nBP >= 2 ~ "Likely benign",
    (PVS1 == 1 & nPS >= 1) |
      (PVS1 == 1 & nPM >= 2) |
      (PVS1 == 1 & nPM == 1 & nPP == 1) |
      (nPS >= 2) |
      (nPS == 1 & nPM >= 3) |
      (nPS == 1 & nPM == 2 & nPP >= 2) |
      (nPS == 1 & nPM == 1 & nPP >= 4) ~ "Pathogenic",
    (PVS1 == 1 & nPM == 1) |
      (PVS1 == 1 & nPP >= 2) |
      (nPS == 1 & nPM >= 1) |
      (nPS == 1 & nPP >= 2) |
      (nPM >= 3) |
      (nPM == 2 & nPP >= 2) |
      (nPM == 1 & nPP >= 4) ~ "Likely pathogenic",
    TRUE ~ "Uncertain significance"
  )

  paste0(
    "InterVar: ", new_class,
    " PVS1=", PVS1,
    " PS=[", paste(PS, collapse = ", "), "]",
    " PM=[", paste(PM, collapse = ", "), "]",
    " PP=[", paste(PP, collapse = ", "), "]",
    " BA1=", BA1,
    " BS=[", paste(BS, collapse = ", "), "]",
    " BP=[", paste(BP, collapse = ", "), "]"
  )
}

# Generate an updated InterVar evidence string for each variant
# incorporating ClinVar-derived PS1 and PM5 evidence
intervar_unique <- intervar_unique %>%
  mutate(
    intervar_updated = pmap_chr(
      list(
        `InterVar: InterVar and Evidence`,
        PS1_new,
        PM5_new
      ),
      update_intervar
    )
  ) %>%
  dplyr::mutate(
    original_class = str_trim(
      str_match(
        `InterVar: InterVar and Evidence`,
        "^InterVar:\\s*(.*?)\\s*PVS1="
      )[, 2]
    ),
    updated_class = str_trim(
      str_match(
        intervar_updated,
        "^InterVar:\\s*(.*?)\\s*PVS1="
      )[, 2]
    ),

    # Compare original and updated InterVar classifications to determine
    # whether the ClinVar evidence changed the final pathogenicity call
    class_changed = original_class != updated_class
  )

# Merge updated InterVar annotations back into the original InterVar
# data frame and replace the existing evidence string when an updated
# version is available
final_intervar_df <- intervar_df %>%
  left_join(intervar_unique %>% dplyr::select(
    `#Chr`, Start, End, Ref, Alt,
    intervar_updated),
    by = join_by(`#Chr`, Start, End, Ref, Alt)) %>%
  dplyr::mutate(`InterVar: InterVar and Evidence` = case_when(
    !is.na(intervar_updated) ~ intervar_updated,
    TRUE ~ `InterVar: InterVar and Evidence`
  ))

# Summarize the number of evidence and classification changes introduced
# by ClinVar-derived PS1/PM5 updates
print(glue::glue("Number of missense variants queried: {nrow(intervar_unique)}"))
print(glue::glue("Number of PS1 updates: {sum(intervar_unique$PS1_old != intervar_unique$PS1_new)}"))

if (sum(intervar_unique$PS1_old != intervar_unique$PS1_new) > 0) {
  intervar_unique %>%
    dplyr::filter(PS1_old != PS1_new) %>%
    dplyr::count(PS1_old, PS1_new)
}

print(glue::glue("Number of PM5 updates: {sum(intervar_unique$PM5_old != intervar_unique$PM5_new)}"))

if (sum(intervar_unique$PM5_old != intervar_unique$PM5_new) > 0) {
  intervar_unique %>%
    dplyr::filter(PM5_old != PM5_new) %>%
    dplyr::count(PM5_old, PM5_new)
}

print(glue::glue("Number of Intervar pathogenicity call updates: {sum(intervar_unique$class_changed == TRUE)}"))

if (sum(intervar_unique$class_changed == TRUE) > 0) {
  intervar_unique %>%
    dplyr::filter(class_changed == TRUE) %>%
    count(original_class, updated_class)
}

# save to output
write_tsv(
  final_intervar_df,
  file.path(results_dir, output_file)
)
