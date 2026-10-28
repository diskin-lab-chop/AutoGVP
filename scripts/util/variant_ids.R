################################################################################
# variant_ids.R
#
# Shared construction of the chrom-pos-ref-alt variant ids used to match the
# VCF, ANNOVAR, InterVar, AutoPVS1 and ClinVar tables.
#
# Chromosome names are written without the "chr" prefix, and the mitochondrial
# contig is written as "M" whether the source calls it chrM, M or MT (ClinVar's
# variant_summary uses MT, UCSC-style VCFs use chrM). This assumes GRCh38, where
# chrM and MT are the same sequence (rCRS); on hg19, UCSC chrM is a different
# sequence and must not be matched to ClinVar's MT.
################################################################################

normalize_chrom <- function(chrom) {
  chrom <- stringr::str_remove(as.character(chrom), "^chr")
  dplyr::if_else(chrom %in% c("M", "MT"), "M", chrom)
}

# Variant id from its parts; positions are written without scientific notation
make_variant_id <- function(chrom, pos, ref, alt) {
  paste(
    normalize_chrom(chrom),
    format(as.numeric(pos), scientific = FALSE, trim = TRUE),
    ref, alt,
    sep = "-"
  )
}

# Normalize the chromosome of an existing chrom-pos-ref-alt id (spaces removed)
normalize_variant_id <- function(id) {
  id <- stringr::str_remove_all(id, " ")
  dplyr::if_else(
    is.na(id),
    NA_character_,
    paste0(normalize_chrom(stringr::str_extract(id, "^[^-]+")), stringr::str_remove(id, "^[^-]+"))
  )
}
