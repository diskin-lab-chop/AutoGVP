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
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "AutoGVP")
input_dir   <- file.path(analysis_dir, "input")
results_dir <- file.path(root_dir, "results")

# create results directory if it does not exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# parse parameters     
option_list <- list(
  make_option(c("--vcf"), type = "character",
              help = "Input filtered and parsed VEP VCF file"),
  make_option(c("--autogvp"), type = "character",
              help = "input AutoGVP annotated file"),
  make_option(c("--output"), type = "character", default = "out",
              help = "output name")) 

opt <- parse_args(OptionParser(option_list = option_list))

# get input files from parameters (read)
input_vcf_file  <-  opt$vcf
input_autogvp_file <- opt$autogvp
output_name   <- opt$output

abridged_out_file <- glue::glue("{output_name}-autogvp-annotated-abridged.tsv")
full_out_file <- glue::glue("{output_name}-autogvp-annotated-full.tsv")


# Read in VEP vcf file
vcf <- read_tsv(file.path(input_dir, input_vcf_file),
                show_col_types = FALSE)

# Remove "[#]" characters from column headers, if present
names(vcf) <- sub(".*]", "", names(vcf))

# create `vcf_id` column to match with AutoGVP output, and rm `ANN` column
vcf <- vcf %>%
  dplyr::mutate(vcf_id = glue::glue("{CHROM}-{POS}-{REF}-{ALT}")) %>%
  dplyr::mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) %>%
  select(-ANN)

# Define names of VEP CSQ fields that will become df column names
csq_cols <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature", 
              "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
              "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM",
              "DISTANCE","STRAND","FLAGS","PICK","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID",
              "CANONICAL","TSL","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM",
              "RefSeq","REFSEQ_MATCH","SOURCE","REFSEQ_OFFSET","GIVEN_REF","USED_REF","BAM_EDIT",
              "GENE_PHENO","SIFT","PolyPhen","DOMAINS","HGVS_OFFSET","HGVSg","AFR_AF","AMR_AF",
              "EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF",
              "gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF",
              "gnomAD_OTH_AF","gnomAD_SAS_AF","CLIN_SIG","SOMATIC","PHENO","PUBMED","CHECK_REF",
              "MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS",
              "CADD_PHRED","CADD_RAW","ALSPAC_AC","ALSPAC_AF","Aloft_pred","BayesDel_noAF_pred",
              "ClinPred_pred","DEOGEN2_pred","Eigen-PC-phred_coding","Eigen-phred_coding",
              "FATHMM_pred","GTEx_V8_gene","GTEx_V8_tissue","Interpro_domain","LIST-S2_pred",
              "LRT_pred","M-CAP_pred","MetaLR_pred","MetaRNN_pred","MetaSVM_pred","MutationAssessor_pred",
              "MutationTaster_pred","PROVEAN_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred",
              "PrimateAI_pred","REVEL_rankscore","REVEL_score","SIFT4G_pred","TWINSUK_AC","TWINSUK_AF",
              "UK10K_AC","UK10K_AF","VEST4_rankscore","VEST4_score","fathmm-MKL_coding_pred",
              "fathmm-XF_coding_pred","gnomAD_exomes_controls_AC","gnomAD_exomes_controls_AF",
              "gnomAD_exomes_controls_AN","gnomAD_exomes_controls_POPMAX_AC","gnomAD_exomes_controls_POPMAX_AF",
              "gnomAD_exomes_controls_POPMAX_AN","gnomAD_exomes_controls_POPMAX_nhomalt",
              "gnomAD_exomes_controls_nhomalt","phastCons100way_vertebrate",
              "phastCons100way_vertebrate_rankscore","phyloP100way_vertebrate",
              "phyloP100way_vertebrate_rankscore","Intervar","Intervar_STATUS")

# Create subset vector of vcf column names to retain
columns_to_retain <- c(names(vcf)[!grepl("CSQ", names(vcf))],
                       "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature", 
                       "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
                       "Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM",
                       "DISTANCE","STRAND","FLAGS","PICK","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID",
                       "CANONICAL","TSL","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM",
                       "RefSeq","REFSEQ_MATCH","SOURCE",
                       "SIFT","PolyPhen","DOMAINS","HGVSg")

# Parse data frame by CSQ subfield and transcript annotation
vcf_separated <- vcf %>%
  
  # Separate `CSQ` transcript annotations into distinct rows (comma separated)
  separate_longer_delim(CSQ, delim = ",") %>%
  
  # separate `CSQ` subfields into unique columns, named in `csq_cols`
  separate_wider_delim(CSQ, "|", names = csq_cols) %>%
  
  dplyr::select(any_of(columns_to_retain))


# Retain one row per variant by prioritizing gene annotation with value "1" in `PICK` column. In KF workflow, criteria for selection include:
# 1) Canonical transcript status
# 2) lowest transcript level support (TSL) category
# 3) Transcript type (protein-coding preferred)
# 4) Highest impact consequence
# 5) CCDC status
# 6) Longest transcript length

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
  dplyr::left_join(vcf[,c("vcf_id", "CSQ")], by = "vcf_id") %>%
  arrange(CHROM, POS)

# Read in autogvp output
autogvp <- read_tsv(file.path(results_dir, input_autogvp_file),
                    show_col_types = FALSE)

# Merge `autogvp` and `vcf_final`
merged_df <- autogvp %>%
  # rm redundant columns from autogvp
  dplyr::select(-CHROM, -ID, -REF, -ALT, -FILTER, -QUAL) %>%
  # join dfs
  dplyr::left_join(vcf_final, by = "vcf_id") %>%
  # rm unecessary columns
  dplyr::select(-any_of(c("var_id", "Otherinfo"))) %>%
  # convert `gnomad_3_1_1_AF_non_cancer` to numeric
  dplyr::mutate(gnomad_3_1_1_AF_non_cancer = case_when(
    gnomad_3_1_1_AF_non_cancer == "." ~ "0",
    TRUE ~ gnomad_3_1_1_AF_non_cancer
  )) %>%
  dplyr::mutate(gnomad_3_1_1_AF_non_cancer = as.numeric(gnomad_3_1_1_AF_non_cancer)) %>%
  # add clinVar link
  dplyr::mutate(ClinVar_link = case_when(
    !is.na(VariationID) ~ paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", VariationID, "/"),
    TRUE ~ NA_character_
  )) %>%
  # rm ensembl transcript/protein IDs from HGVSc/p columns
  dplyr::mutate(HGVSc = str_split(HGVSc, ":", simplify = T)[,2],
                HGVSp = str_split(HGVSp, ":", simplify = T)[,2]) %>%
  arrange(CHROM, START)

# Create list of columns to include in abridged output
abridged_cols <- c("CHROM", "START", "REF", "ALT", "ID",
                   "SYMBOL", "Consequence", "HGVSg", "HGVSc", "HGVSp", 
                   "final_call", "Reasoning_for_call", "Stars",
                   "ClinVar_ClinicalSignificance",
                   "Intervar_evidence", "gnomad_3_1_1_AF_non_cancer")
  
# Generate abridged output
abridged_df <- merged_df %>%
  dplyr::select(all_of(abridged_cols)) %>%
  write_tsv(file.path(results_dir, abridged_out_file))

# Generate comprehensive output
full_df <- merged_df %>%
  dplyr::relocate(all_of(abridged_cols)) %>%
  write_tsv(file.path(results_dir, full_out_file))
