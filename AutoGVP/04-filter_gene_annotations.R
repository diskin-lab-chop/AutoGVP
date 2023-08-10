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

# set up directories
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
    help = "Input filtered and parsed VEP VCF file"
  ),
  make_option(c("--autogvp"),
    type = "character",
    help = "input AutoGVP annotated file"
  ),
  make_option(c("--output"),
    type = "character", default = "out",
    help = "output name"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# get input files from parameters
input_vcf_file <- opt$vcf
input_autogvp_file <- opt$autogvp
output_name <- opt$output

# Define output file variables
abridged_out_file <- glue::glue("{output_name}-autogvp-annotated-abridged.tsv")
full_out_file <- glue::glue("{output_name}-autogvp-annotated-full.tsv")

# Define path for VCF file CSQ field names
csq_fields <- file.path(results_dir, glue::glue("{output_name}.filtered_csq_subfields.tsv"))

# Read in VEP vcf file
vcf <- read_tsv(input_vcf_file,
  show_col_types = FALSE
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
  dplyr::select(any_of(columns_to_retain))


# Retain one row per variant by prioritizing gene annotation with value "1" in `PICK` column. In KF workflow, criteria for selection include:
# 1) Canonical transcript status
# 2) lowest transcript level support (TSL) category
# 3) Transcript type (protein-coding preferred)
# 4) Highest impact consequence
# 5) CCDC status
# 6) Longest transcript length

# Retain only transcripts where PICK == 1
vcf_pick <- vcf_separated %>%
  dplyr::filter(PICK == "1")

# In some instances, entrezGene annotation is chosen when Ensembl has a comparable annotation. We can manually replace entrezGene with ensembl annotation in these cases
# Extract ensembl annotations:
vcf_pick_ensembl <- vcf_pick %>%
  dplyr::filter(grepl("ENS", Feature))

# Identify pick annotations for which a canonical ensembl annotation is also available:
vcf_pick_other <- vcf_pick %>%
  dplyr::filter(!grepl("ENS", Feature)) %>%
  dplyr::mutate(vcf_id_gene_csq = glue::glue("{vcf_id}-{SYMBOL}-{Consequence}"))

# Identify pick annotations for which an ensembl annotation to the same gene and consequence is available:
vcf_pick_other_ensembl <- vcf_separated %>%
  dplyr::filter(grepl("ENS", Feature)) %>%
  dplyr::filter(glue::glue("{vcf_id}-{SYMBOL}-{Consequence}") %in% vcf_pick_other$vcf_id_gene_csq) %>%
  group_by(vcf_id) %>%
  arrange(CANONICAL) %>%
  dplyr::slice_tail(n = 1)

# Bind df rows to create `vcf_final` df:
vcf_final <- vcf_pick_other %>%
  dplyr::select(-vcf_id_gene_csq) %>%
  dplyr::filter(!vcf_id %in% vcf_pick_other_ensembl$vcf_id) %>%
  bind_rows(vcf_pick_other_ensembl, vcf_pick_ensembl) %>%
  # Add `csq` column with all annotation info
  dplyr::left_join(vcf[, c("vcf_id", "CSQ")], by = "vcf_id") %>%
  arrange(CHROM, POS)


# Read in autogvp output
autogvp <- read_tsv(input_autogvp_file,
  show_col_types = FALSE
)

# Parse Sample column, if present
if ("Sample" %in% names(autogvp)) {
  autogvp <- autogvp %>%
    tidyr::separate_wider_delim(Sample, delim = ":", names = c("GT", "AD", "DP", "GQ"), too_many = "drop")
}

# Merge `autogvp` and `vcf_final`
merged_df <- autogvp %>%
  # rm redundant columns from autogvp
  dplyr::select(-any_of(c(
    "CHROM", "POS", "ID", "REF",
    "ALT", "FILTER", "QUAL",
    "Interpro_domain"
  ))) %>%
  # join dfs
  dplyr::left_join(vcf_final, by = "vcf_id") %>%
  # rm unecessary columns
  dplyr::select(-any_of(c("var_id", "Otherinfo"))) %>%
  # add clinVar link
  dplyr::mutate(ClinVar_link = case_when(
    !is.na(VariationID) ~ paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", VariationID, "/"),
    TRUE ~ NA_character_
  )) %>%
  # coordinate-sort
  arrange(CHROM, POS)

# Change `gnomad_3_1_1_AF_non_cancer` to numeric, when present
if ("gnomad_3_1_1_AF_non_cancer" %in% names(merged_df)) {
  merged_df <- merged_df %>%
    # convert `gnomad_3_1_1_AF_non_cancer` to numeric
    dplyr::mutate(gnomad_3_1_1_AF_non_cancer = case_when(
      gnomad_3_1_1_AF_non_cancer == "." ~ "0",
      TRUE ~ gnomad_3_1_1_AF_non_cancer
    )) %>%
    dplyr::mutate(gnomad_3_1_1_AF_non_cancer = as.numeric(gnomad_3_1_1_AF_non_cancer))
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
    )
}

# Define `coacross()` function to coalesce across multiple df columns
coacross <- function(...) {
  coalesce(!!!across(...))
}

# Coalesce rsIDs across multiple columns, when present
id_df <- merged_df %>%
  dplyr::select(any_of(c("ID", "avsnp147", "Existing_variation"))) %>%
  dplyr::mutate_all(funs(ifelse(grepl("rs", .), ., NA_character_))) %>%
  dplyr::mutate(rsID = coacross())

# Add coalesced rsID to merged_df, and remove other ID columns
merged_df <- merged_df %>%
  dplyr::mutate(rsID = id_df$rsID) %>%
  select(-any_of(c("ID", "avsnp147", "Existing_variation")))

# read in output file column names tsv
colnames <- read_tsv(file.path(input_dir, "output_colnames.tsv"))

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
  write_tsv(file.path(results_dir, abridged_out_file))

# Generate comprehensive output
full_df <- merged_df %>%
  dplyr::relocate(all_of(abridged_cols)) %>%
  write_tsv(file.path(results_dir, full_out_file))
