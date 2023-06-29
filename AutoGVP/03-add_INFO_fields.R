################################################################################
# 03-add_INFO_fields.R
# written by Ammar Naqvi
#
# This script annotates INFO column variables with variants that have been 
# processed and re-called by AutoGVP
#
# usage: Rscript  03-add_INFO_fields.R  --input_table <tsv file> 
#                                       --vcf_file <vcf file> 
#                                       --out_prefix <string>
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
  library("vcfR")
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "AutoGVP")
results_dir <- file.path(root_dir, "results")
input_dir   <- file.path(analysis_dir, "input")

# parse parameters     
option_list <- list(
                    make_option(c("--input_table"), type = "character",
                                help = "Output tab from AutoGVP"),
                    make_option(c("--vcf_file"), type = "character",
                                help = "vcf file used"),
                    make_option(c("--out_prefix"), type = "character",
                                help = "output prefix for result file ")
              )

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (reqd)
input_tab  <- opt$input_table
vcf_file   <- opt$vcf_file
output_prefix <- opt$out_prefix

output_tab_file  <- file.path(results_dir, paste0(output_prefix,"_annotations_report.full.tsv"))

input_tab <- vroom(input_tab, show_col_types = TRUE)
vcfR_df <- read.vcfR(vcf_file, verbose = FALSE )

## get INFO fields and convert to df for each variable
INFO_df <- INFO2df(vcfR_df)

## combine new columns with abridged/summary table and write to file
vcf_added_columns_df <- cbind(input_tab,INFO_df) %>% 
                        dplyr::rename("Hugo_Symbol"=Ref.Gene) %>% 
                        write_tsv(output_tab_file)
