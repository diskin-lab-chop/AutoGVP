################################################################################
# update_intervar.R
# written by Ryan Corbett
#
# This script loads the specified intervar df and updates PS1 and PM5 criteria
# based on clinvar entries in supplied resolved clinvar
# interpretations
#
# usage: update_intervar.R --intervar_file <intervar file>
#                          --clinvar_file <resolved clinvar file>
#                          --clinvar_hgvs4_file <clinvar hgvs4 variation file>
#                          --vcf_file <VEP-annotated vcf>
#
################################################################################

options(scipen = 999)

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
  make_option(c("--vcf_file"),
    type = "character",
    help = "VEP-annotated vcf (RefSeq transcripts used for PS1/PM5 matching)"
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
vcf_file <- opt$vcf_file
results_dir <- opt$outdir
output_file <- glue::glue("{basename(intervar_file)}_updated")

# Load InterVar output and ensure chromosome values are stored as character
# strings and start is numeric to avoid type mismatches during downstream joins
intervar_df <- read_tsv(intervar_file, show_col_types = FALSE) %>%
  dplyr::mutate(
    `#Chr` = as.character(`#Chr`),
    Start = as.numeric(Start)
  )

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

# Read the "##" meta-information lines of a (optionally gzipped) vcf
read_vcf_header <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  header <- character()
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0 || !startsWith(line, "##")) break
    header <- c(header, line)
  }
  header
}

# VEP CSQ subfield names, from the CSQ INFO header line ("Format: Allele|Consequence|...")
csq_header <- str_subset(read_vcf_header(vcf_file), "^##INFO=<ID=CSQ,")
if (length(csq_header) == 0) {
  stop(glue::glue("No VEP CSQ header found in {vcf_file}; a VEP-annotated vcf is required"))
}
csq_cols <- str_split_1(str_match(csq_header[1], "Format: ([^\"]+)")[, 2], "\\|")
missing_csq_cols <- setdiff(c("Allele", "Consequence", "SYMBOL", "Feature", "HGVSp"), csq_cols)
if (length(missing_csq_cols) > 0) {
  stop(glue::glue("VEP CSQ field(s) missing from {vcf_file}: {paste(missing_csq_cols, collapse = ', ')}"))
}

# Protein change of each variant on every RefSeq (NM_) transcript where VEP calls it a
# missense, so that ClinVar records can be compared on the same transcript.
# Multi-allelic ALTs are split and each CSQ entry is kept for its own allele.
vep_missense_df <- vroom::vroom(
  vcf_file,
  delim = "\t",
  comment = "##",
  col_select = c(`#CHROM`, POS, REF, ALT, INFO),
  col_types = vroom::cols(.default = "c"),
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    chr_norm = str_remove(`#CHROM`, "^chr"),
    Start = as.numeric(POS),
    Ref = REF,
    Alt = ALT,
    CSQ = str_match(INFO, "(?:^|;)CSQ=([^;]+)")[, 2]
  ) %>%
  dplyr::filter(!is.na(CSQ)) %>%
  tidyr::separate_longer_delim(Alt, delim = ",") %>%
  tidyr::separate_longer_delim(CSQ, delim = ",") %>%
  tidyr::separate_wider_delim(CSQ, delim = "|", names = csq_cols) %>%
  dplyr::filter(
    Allele == Alt,
    str_detect(Consequence, "missense_variant"),
    str_detect(Feature, "^NM_"),
    HGVSp != ""
  ) %>%
  dplyr::mutate(
    hgvsp_parts = str_match(
      convert_hgvsp_to_one_letter(str_remove(HGVSp, "^[^:]*:")),
      "^p\\.([A-Z])(\\d+)([A-Z])$"
    ),
    ref_aa = hgvsp_parts[, 2],
    aa_pos = as.integer(hgvsp_parts[, 3]),
    alt_aa = hgvsp_parts[, 4],
    gene = SYMBOL,
    transcript = Feature,
    transcript_acc = str_remove(Feature, "\\.\\d+$")
  ) %>%
  dplyr::filter(!is.na(aa_pos)) %>%
  dplyr::distinct(chr_norm, Start, Ref, Alt, gene, transcript, transcript_acc, ref_aa, aa_pos, alt_aa)

print(glue::glue(
  "VEP RefSeq missense annotations: {nrow(vep_missense_df)} across ",
  "{nrow(dplyr::distinct(vep_missense_df, chr_norm, Start, Ref, Alt))} variants"
))

# PS1/PM5 are evaluated at the amino acid (codon) level: a P/LP ClinVar missense at the
# same residue of the same RefSeq transcript supports the criteria, regardless of which
# nucleotide of the codon is changed.
#   PS1: same amino acid change as the InterVar variant, caused by a different nucleotide change
#   PM5: different missense amino acid change at the same residue
# Residues are compared on a shared transcript, because residue numbers differ between
# isoforms; codons are defined by transcript coordinates, so codons split across exons are
# handled. Transcripts are matched by accession regardless of version, and the reference
# residue (e.g. G in p.G172S) must also match, to guard against numbering changes between
# transcript versions.

# subset intervar df for missense variants
intervar_missense_df <- intervar_df %>%
  dplyr::filter(ExonicFunc.refGene == "nonsynonymous SNV") %>%
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
    chr_norm = str_remove(`#Chr`, "^chr"),
    var_id = paste(chr_norm, format(Start, scientific = FALSE, trim = TRUE), Ref, Alt, sep = "-")
  ) %>%
  dplyr::select(
    `#Chr`, Start, End, Ref, Alt, chr_norm, var_id, Ref.Gene,
    `clinvar: Clinvar`,
    `InterVar: InterVar and Evidence`,
    PS1_old, PM5_old
  ) %>%
  dplyr::rename(ClinicalSignificance_old_clinvar = `clinvar: Clinvar`)

# protein change of each missense variant on every RefSeq transcript annotated by VEP;
# variants without a VEP RefSeq missense annotation cannot be evaluated
intervar_missense_aa_df <- intervar_missense_df %>%
  dplyr::distinct(`#Chr`, Start, End, Ref, Alt, chr_norm, var_id) %>%
  dplyr::inner_join(
    vep_missense_df,
    by = c("chr_norm", "Start", "Ref", "Alt"),
    relationship = "one-to-many"
  )

# P/LP records in the resolved ClinVar interpretations
# (a word boundary is used so "Conflicting classifications of pathogenicity" is not
# treated as a P/LP record)
clinvar_plp_df <- clinvar_df %>%
  dplyr::filter(
    str_detect(ClinicalSignificance, regex("\\bpathogenic\\b", ignore_case = TRUE)),
    !str_detect(ClinicalSignificance, regex("benign", ignore_case = TRUE))
  ) %>%
  dplyr::transmute(
    ClinVar_VariationID = VariationID,
    clinvar_chr = str_remove(chr_clinvar, "^chr"),
    clinvar_pos = Start_clinvar,
    clinvar_var_id = paste(
      clinvar_chr,
      format(Start_clinvar, scientific = FALSE, trim = TRUE),
      Ref_clinvar, Alt_clinvar,
      sep = "-"
    )
  ) %>%
  dplyr::distinct()

# Load ClinVar HGVS4Variation file, which provides protein-level HGVS
# annotations (ProteinChange) linked to ClinVar VariationIDs.
# This file is typically tens of GB uncompressed, but only rows for P/LP VariationIDs
# on a transcript of an InterVar missense variant are needed.
# Pre-filter with awk (a fast, single C-level pass over the decompressed
# stream) so that only the small matching subset ever reaches vroom,
# instead of loading and indexing the entire file in R.
plp_ids <- unique(clinvar_plp_df$ClinVar_VariationID)
query_transcripts <- unique(intervar_missense_aa_df$transcript_acc)

# If there is nothing to look up in the HGVS4Variation file, no PS1/PM5 evidence
# can be updated. Write the input through unchanged and exit early rather
# than handing vroom an empty pipe further down.
if (length(plp_ids) == 0 || length(query_transcripts) == 0) {
  print("No missense variants or P/LP ClinVar records to evaluate; skipping PS1/PM5 updates.")
  write_tsv(intervar_df, file.path(results_dir, output_file))
  quit(save = "no", status = 0)
}

decompress_cmd <- "gzip -cd"

ids_file <- tempfile()
writeLines(sprintf("%.0f", plp_ids), ids_file)
transcripts_file <- tempfile()
writeLines(query_transcripts, transcripts_file)
# query transcripts found on any ClinVar record (P/LP or not), written by the awk filter
seen_transcripts_file <- tempfile()

hgvs4_cols <- c(
  "Symbol", "GeneID", "VariationID", "AlleleID", "Type", "Assembly",
  "NucleotideExpression", "NucleotideChange", "ProteinExpression",
  "ProteinChange", "UsedForNaming", "Submitted", "OnRefSeqGene"
)

# VariationID is column 3 and NucleotideExpression (e.g. NM_000059.4:c.10153C>T) is
# column 7 of the HGVS4Variation file; keep only data lines (drop the leading "#..."
# comment/header lines) whose VariationID is a P/LP record and whose transcript accession,
# without version, is a transcript of an InterVar missense variant. Query transcripts that
# appear on any record are written to `seen_transcripts_file`.
filter_cmd <- sprintf(
  "%s %s | awk -F'\\t' -v seen_file=%s 'FILENAME == ARGV[1] { gsub(/\\r$/, \"\"); ids[$1]; next } FILENAME == ARGV[2] { gsub(/\\r$/, \"\"); tx[$1]; next } { gsub(/\\r$/, \"\") } !/^#/ { acc = $7; sub(/[.:].*$/, \"\", acc); if (!(acc in tx)) next; if (!(acc in seen)) { seen[acc]; print acc > seen_file } if ($3 in ids) print }' %s %s -",
  decompress_cmd,
  shQuote(clinvar_hgvs4_file),
  shQuote(seen_transcripts_file),
  shQuote(ids_file),
  shQuote(transcripts_file)
)

hgvs4_variation_df <- vroom::vroom(
  pipe(filter_cmd),
  delim = "\t",
  col_names = hgvs4_cols,
  col_select = c(VariationID, Assembly, NucleotideExpression, ProteinChange),
  col_types = c(VariationID = "n", Assembly = "c", NucleotideExpression = "c", ProteinChange = "c"),
  show_col_types = FALSE
)

clinvar_transcripts <- if (file.exists(seen_transcripts_file)) readLines(seen_transcripts_file) else character()
unlink(c(ids_file, transcripts_file, seen_transcripts_file))

# Variants that cannot be compared with ClinVar on any transcript keep InterVar's PS1/PM5
# evidence: those without a VEP RefSeq missense annotation, and those none of whose
# transcripts appear in ClinVar (e.g. a transcript ClinVar has replaced since the VEP
# cache was built)
no_vep_df <- intervar_missense_df %>%
  dplyr::anti_join(intervar_missense_aa_df, by = "var_id")
if (nrow(no_vep_df) > 0) {
  print(glue::glue(
    "{dplyr::n_distinct(no_vep_df$var_id)} missense variant(s) without a VEP RefSeq missense ",
    "annotation; InterVar PS1/PM5 kept: {paste(sort(unique(no_vep_df$Ref.Gene)), collapse = ', ')}"
  ))
}

no_clinvar_tx_df <- intervar_missense_aa_df %>%
  dplyr::group_by(var_id) %>%
  dplyr::filter(!any(transcript_acc %in% clinvar_transcripts)) %>%
  dplyr::ungroup()
if (nrow(no_clinvar_tx_df) > 0) {
  print(glue::glue(
    "{dplyr::n_distinct(no_clinvar_tx_df$var_id)} missense variant(s) with no transcript in ClinVar; ",
    "InterVar PS1/PM5 kept: {paste(sort(unique(no_clinvar_tx_df$gene)), collapse = ', ')}"
  ))
}

# Retain only protein-level annotations from the HGVS4Variation file.
# Assembly == "na" corresponds to protein annotations rather than
# genomic or transcript-level representations.
# Convert amino acid abbreviations to match InterVar formatting
# (e.g. p.Trp507Arg -> p.W507R) and retain missense substitutions only; this
# excludes nonsense (p.E285*), synonymous (p.G10=), start-loss/unknown (p.M1?),
# frameshift, extension and in-frame indel consequences.
clinvar_protein_df <- hgvs4_variation_df %>%
  dplyr::filter(
    Assembly == "na",
    ProteinChange != "-"
  ) %>%
  dplyr::mutate(
    transcript_clinvar = str_remove(NucleotideExpression, ":.*$"),
    transcript_acc = str_remove(transcript_clinvar, "\\.\\d+$"),
    ProteinChange = convert_hgvsp_to_one_letter(ProteinChange),
    clinvar_parts = str_match(ProteinChange, "^p\\.([A-Z])(\\d+)([A-Z])$"),
    clinvar_ref_aa = clinvar_parts[, 2],
    aa_pos = as.integer(clinvar_parts[, 3]),
    clinvar_alt_aa = clinvar_parts[, 4]
  ) %>%
  dplyr::filter(!is.na(aa_pos)) %>%
  dplyr::distinct(VariationID, transcript_clinvar, transcript_acc, clinvar_ref_aa, aa_pos, clinvar_alt_aa)

# Match each InterVar missense annotation to P/LP ClinVar missense records at the same
# residue of the same transcript (any version), then evaluate PS1 and PM5 for each match.
intervar_matches_df <- intervar_missense_aa_df %>%
  dplyr::inner_join(
    clinvar_protein_df,
    by = c("transcript_acc", "aa_pos", "ref_aa" = "clinvar_ref_aa"),
    relationship = "many-to-many"
  ) %>%
  dplyr::inner_join(
    clinvar_plp_df,
    by = c("VariationID" = "ClinVar_VariationID"),
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    # a record for the same nucleotide change is the variant itself
    different_nt = var_id != clinvar_var_id,
    PS1_support = different_nt & alt_aa == clinvar_alt_aa,
    PM5_support = different_nt & alt_aa != clinvar_alt_aa
  )

# Report supporting matches for which ClinVar does not list the record on VEP's transcript
# version, i.e. the residue was compared on another version of the transcript
version_mismatch_df <- intervar_matches_df %>%
  dplyr::filter(PS1_support | PM5_support) %>%
  dplyr::group_by(var_id, VariationID, transcript_acc) %>%
  dplyr::filter(!any(transcript == transcript_clinvar)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(var_id, VariationID, transcript, transcript_clinvar)

if (nrow(version_mismatch_df) > 0) {
  print(glue::glue(
    "{nrow(dplyr::distinct(version_mismatch_df, var_id, VariationID))} PS1/PM5 match(es) ",
    "made on a different transcript version than VEP's:"
  ))
  print(version_mismatch_df)
}

# Aggregate evidence across all ClinVar matches associated with the
# same InterVar variant (any transcript). A single supporting ClinVar record
# is sufficient to activate PS1 or PM5 for the variant.
variant_summary <- intervar_matches_df %>%
  dplyr::group_by(`#Chr`, Start, End, Ref, Alt) %>%
  dplyr::summarise(
    PS1_new = as.integer(any(PS1_support)),
    PM5_new = as.integer(any(PM5_support)),
    PS1_ClinVarIDs = paste(unique(VariationID[PS1_support]), collapse = ";"),
    PM5_ClinVarIDs = paste(unique(VariationID[PM5_support]), collapse = ";"),
    .groups = "drop"
  )

# One row per InterVar variant (InterVar reports a single evidence string per variant,
# although a variant can have several transcript annotations), with the aggregated
# PS1/PM5 evidence assignments calculated from ClinVar
intervar_unique <- intervar_missense_df %>%
  dplyr::distinct(
    `#Chr`, Start, End, Ref, Alt,
    ClinicalSignificance_old_clinvar,
    `InterVar: InterVar and Evidence`,
    PS1_old, PM5_old
  ) %>%
  dplyr::left_join(
    variant_summary,
    by = c("#Chr", "Start", "End", "Ref", "Alt")
  ) %>%
  # variants without a P/LP record at the residue have no new evidence
  dplyr::mutate(
    PS1_new = coalesce(PS1_new, 0L),
    PM5_new = coalesce(PM5_new, 0L)
  )

# function to update the PS1 and PM5 values in the `InterVar: InterVar and Evidence`
# string. Only the PS and PM evidence vectors are modified; the InterVar
# classification is NOT recalculated here. Final classification (including
# exclusion of PP5/BP6 and handling of conflicting evidence) is performed
# in 02-annotate_variants.R from the evidence vectors, so the leading class
# label in the updated string should not be interpreted as the final call.
update_intervar <- function(intervar_string, PS1_new, PM5_new) {
  PS <- str_match(intervar_string, "PS=\\[([^]]+)\\]")[, 2] |>
    str_split(",\\s*") |>
    unlist() |>
    as.numeric()

  PM <- str_match(intervar_string, "PM=\\[([^]]+)\\]")[, 2] |>
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

  intervar_string %>%
    str_replace("PS=\\[[^]]+\\]", paste0("PS=[", paste(PS, collapse = ", "), "]")) %>%
    str_replace("PM=\\[[^]]+\\]", paste0("PM=[", paste(PM, collapse = ", "), "]"))
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
  )

# Merge updated InterVar annotations back into the original InterVar
# data frame and replace the existing evidence string when an updated
# version is available
final_intervar_df <- intervar_df %>%
  left_join(
    intervar_unique %>% dplyr::select(
      `#Chr`, Start, End, Ref, Alt,
      intervar_updated
    ),
    by = join_by(`#Chr`, Start, End, Ref, Alt)
  ) %>%
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

# save to output
write_tsv(
  final_intervar_df,
  file.path(results_dir, output_file)
)
