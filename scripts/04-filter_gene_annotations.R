################################################################################
# 04-parse_gene_annotation_merge.R
# written by Ryan Corbett
#
# This scripts reads in filtered, subfield-parsed VEP VCF file, and further parses
# VEP `CSQ` subfield such that there is a unique gene annotation row for each
# variant. This data frame is then merged with AutoGVP output, and full
# and abbreviated results files are written to output.
#
# usage: Rscript 04-parse_gene_annotation_merge.R --vcf <vcf file>
#                                       --autogvp <autogvp file>
#                                       --output <output prefix>
#
################################################################################

# Load required packages
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
    help = "Input filtered and parsed VEP VCF file"
  ),
  make_option(c("--autogvp"),
    type = "character",
    help = "input AutoGVP annotated file"
  ),
  make_option(c("--output"),
    type = "character", default = "out",
    help = "output name"
  ),
  make_option(c("--outdir"),
    type = "character", default = "results",
    help = "output directory"
  ),
  make_option(c("--default_colnames"),
    type = "character", default = "data/output_colnames_default.tsv",
    help = "default output colnames"
  ),
  make_option(c("--custom_colnames"),
              type = "character", default = NULL,
              help = "user-defined output colnames"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# get input files from parameters
input_vcf_file <- opt$vcf
input_autogvp_file <- opt$autogvp
output_name <- opt$output
results_dir <- opt$outdir
default_colnames_file <- opt$default_colnames
custom_colnames_file <- opt$custom_colnames

# create results directory if it does not exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define output file variables
abridged_out_file <- glue::glue("{output_name}-autogvp-annotated-abridged.tsv")
full_out_file <- glue::glue("{output_name}-autogvp-annotated-full.tsv")

# Define path for VCF file CSQ field names
csq_fields <- file.path(results_dir, glue::glue("{output_name}.filtered_csq_subfields.tsv"))

# Read in VEP vcf file
vcf <- read_tsv(input_vcf_file,
  show_col_types = FALSE,
  guess_max = 10000
)

# Remove "[#]" characters from column headers, if present
names(vcf) <- sub(".*]", "", names(vcf))

# create `vcf_id` column to match with AutoGVP output, and rm `ANN` column
vcf <- vcf %>%
  dplyr::mutate(vcf_id = glue::glue("{CHROM}-{POS}-{REF}-{ALT}")) %>%
  dplyr::mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) %>%
  select(-any_of(c("ANN")))

# Read in VEP CSQ field names
csq_cols <- read_lines(csq_fields)

# Create subset vector of vcf column names to retain
columns_to_retain <- c(
  names(vcf)[!grepl("CSQ|DP", names(vcf))],
  "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
  "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
  "Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM",
  "DISTANCE", "STRAND", "FLAGS", "PICK", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID",
  "CANONICAL", "TSL", "CCDS", "ENSP", "Interpro_domain", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM",
  "RefSeq", "REFSEQ_MATCH", "SOURCE",
  "SIFT", "PolyPhen", "DOMAINS", "HGVSg"
)

# Parse data frame by CSQ field and transcript annotation
vcf_separated <- vcf %>%
  # Separate `CSQ` transcript annotations into distinct rows (comma separated)
  separate_longer_delim(CSQ, delim = ",") %>%
  # separate `CSQ` fields into unique columns, named in `csq_cols`
  separate_wider_delim(CSQ, "|", names = csq_cols) %>%
  # Select columns to retain
  dplyr::select(any_of(columns_to_retain)) %>%
  dplyr::select(-any_of(c(
    "CHROM", "POS", "ID", "REF",
    "ALT", "FILTER", "QUAL",
    "Interpro_domain"
  )))

# Read in autogvp output
autogvp <- read_tsv(input_autogvp_file,
  show_col_types = FALSE,
  guess_max = 10000
)

# Parse Sample column, if present
if ("Sample" %in% names(autogvp)) {
  
  # write function to parse SAMPLE format fields based on 
  parse_sample <- function(fmt, smp) {
    ff <- strsplit(fmt, ":")[[1]]
    sv <- strsplit(smp, ":")[[1]]
    res <- setNames(as.list(sv), ff)
    # pad missing fields with NA
    missing <- setdiff(ff, names(res))
    res[missing] <- NA
    as.data.frame(res)
  }
  
  autogvp_expanded <- dplyr::bind_rows(
    purrr::map2(autogvp$FORMAT, autogvp$Sample, parse_sample)
  )
  
  autogvp <- bind_cols(autogvp, autogvp_expanded) %>%
    tidyr::separate_wider_delim(AD, delim = ",", names = c("AD_ref", "AD_alt"), too_many = "drop")
  
}

# Merge `autogvp` and `vcf_final`
merged_df <- autogvp %>%
  # join dfs
  dplyr::left_join(vcf_separated, by = c("vcf_id", "Feature")) %>%
  dplyr::left_join(vcf[, c("vcf_id", "CSQ")], by = "vcf_id") %>%
  # rm unecessary columns
  dplyr::select(-any_of(c("var_id", "Otherinfo"))) %>%
  # add clinVar link
  dplyr::mutate(ClinVar_link = case_when(
    !is.na(VariationID) ~ paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", VariationID, "/"),
    TRUE ~ NA_character_
  ))

# Change `gnomad_3_1_1_AF_non_cancer_popmax` to numeric, when present
if ("gnomad_3_1_1_AF_non_cancer_popmax" %in% names(merged_df)) {
  merged_df <- merged_df %>%
    # convert `gnomad_3_1_1_AF_non_cancer_popmax` to numeric
    dplyr::mutate(gnomad_3_1_1_AF_non_cancer_popmax = case_when(
      gnomad_3_1_1_AF_non_cancer_popmax == "." ~ "0",
      TRUE ~ gnomad_3_1_1_AF_non_cancer_popmax
    )) %>%
    dplyr::mutate(gnomad_3_1_1_AF_non_cancer_popmax = as.numeric(gnomad_3_1_1_AF_non_cancer_popmax))
}

# Reformat HGVSc and HGVSp columns, when present, to remove gene IDs
if ("HGVSc" %in% names(merged_df)) {
  merged_df <- merged_df %>%
    # rm ensembl transcript/protein IDs from HGVSc columns
    dplyr::mutate(
      HGVSc = str_split(HGVSc, ":", simplify = T)[, 2],
    )
}

if ("HGVSp" %in% names(merged_df)) {
  merged_df <- merged_df %>%
    # rm ensembl transcript/protein IDs from HGVSp columns
    dplyr::mutate(
      HGVSp = str_split(HGVSp, ":", simplify = T)[, 2],
    ) %>%
    # replace `%3D` symbol with `=` in synonymous variant HGVSp
    dplyr::mutate(
      HGVSp = str_replace(HGVSp, "%3D", "=")
    )
}

split_and_unique <- function(string) {
  # Split the string by both ";" and ","
  split_values <- unlist(strsplit(string, "[;,]"))
  # Return only unique values
  unique_values <- unique(split_values[!is.na(split_values)])
  unique_values <- paste(unique_values, collapse = ";")
  unique_values <- unlist(unique_values)

  return(unique_values)
}

# Coalesce rsIDs across multiple columns, when present
id_df <- merged_df %>%
  dplyr::select(any_of(c("ID", "avsnp147", "Existing_variation"))) %>%
  dplyr::mutate(across(everything(), ~ ifelse(.x %in% c(".", NA_character_), "", .x))) %>%
  rowwise() %>%
  mutate(variantID = glue::glue(paste(c_across(everything()), collapse = ";"))) %>%
  ungroup() %>%
  dplyr::mutate(variantID = unlist(lapply(variantID, split_and_unique))) %>%
  dplyr::mutate(variantID = case_when(
    variantID == "NA" ~ "",
    TRUE ~ variantID
  )) %>%
  dplyr::mutate(variantID = str_replace_all(variantID, "~", "")) %>%
  dplyr::mutate(variantID = str_replace(variantID, "^;", ""))

# Add coalesced rsID to merged_df, and remove other ID columns
merged_df <- merged_df %>%
  dplyr::mutate(variantID = id_df$variantID) %>%
  select(-any_of(c("ID", "avsnp147", "Existing_variation")))

# read in default output column names tsv
default_colnames <- read_tsv(default_colnames_file,
                             show_col_types = FALSE
)

# if custom output colnames provided, append to default colnames
if (!is.null(custom_colnames_file)) {
  
  custom_colnames <- read_tsv(custom_colnames_file,
                              show_col_types = FALSE)
  
  if (length(names(custom_colnames)) != 3 & all(names(custom_colnames) != c("Column_name", "Rename", "Abridged"))){
    
    stop("Error: custom_colnames should contain three columns with names 'Column_name', 'Rename', 'Abridged')")
    
  }
  
  colnames <- default_colnames %>%
    # remove col_name from default if also in custom to ensure abridged status is defined by user 
    dplyr::filter(!Column_name %in% custom_colnames$Column_name) %>%
    bind_rows(custom_colnames)

} else {
  
  colnames <- default_colnames
  
}

# Subset and reorder output columns based on inclusion and order in `colnames`
merged_df <- merged_df %>%
  dplyr::select(any_of(colnames$Column_name))

# filter `colnames` for columns present in `merged_df`
colnames <- colnames %>%
  dplyr::filter(Column_name %in% names(merged_df))

# rename columns according to `colnames` Rename column
merged_df <- merged_df %>%
  rename_at(vars(colnames$Column_name), ~ colnames$Rename)

# Create list of columns to include in abridged output
abridged_cols <- colnames %>%
  filter(Abridged == T) %>%
  pull(Rename)

# Generate abridged output
abridged_df <- merged_df %>%
  dplyr::select(any_of(abridged_cols)) %>%
  dplyr::arrange(chr, start) %>%
  write_tsv(file.path(results_dir, abridged_out_file))

# Generate comprehensive output
full_df <- merged_df %>%
  dplyr::relocate(all_of(abridged_cols)) %>%
  dplyr::arrange(chr, start) %>%
  write_tsv(file.path(results_dir, full_out_file))
