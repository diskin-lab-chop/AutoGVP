################################################################################
# 01-annotate_variants.R
# written by Ammar Naqvi
#
# This script annotates variants based on clinVar, and recomputes score and call 
# for intervar and autoPVS1
#
# usage: Rscript 01-annotate_variants.R --vcf <vcf file> 
#                                       --intervar <intervar file> 
#                                       --autopvs1 <autopvs1 file>
#                                       --clinvar  'yyyymmdd'
#                                       --gnomad_variable 'gnomad_3_1_1_AF_non_cancer'
#                                       --gnomad_af <numeric>
#                                       --variant_depth <integer>
#                                       --variant_af <numeric>
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  library("vroom")
} )

#Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses")
input_dir   <- file.path(analysis_dir, "input")

# parse parameters     
option_list <- list(
  make_option(c("--vcf"), type = "character",
              help = "Input vcf file with clinVar annotations"),
  make_option(c("--intervar"), type = "character",
              help = "input intervar file"),
  make_option(c("--clinvar"), type = "character",
              help = "specific clinVar file (format: clinvar_20211225.vcf.gz)"), ## https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022clinvar_20211225.vcf.gz 
  
  make_option(c("--gnomad_variable"), type = "character",default = "gnomAD_genome_ALL",
              help = "gnomAD variable"),
  make_option(c("--gnomad_af"), type = "numeric", default = 0.001,
              help = "genomAD AF filter"),
  make_option(c("--variant_depth"), type = "integer",default = 15,
              help = "variant depth filter"),
  make_option(c("--variant_af"), type = "numeric", default = .2,
              help = "variant AF cut-off")) 

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (reqd)
input_vcf_file <- opt$vcf
input_clinVar_file  <-  opt$clinvar
input_intervar_file <- opt$intervar

## for testing purposes 
input_clinVar_file  <- "/Users/naqvia/d3b_codes/pathogenicity-assessment/analyses/input/clinvar_20220403.vcf.gz" 
input_intervar_file <- "/Users/naqvia/d3b_codes/pathogenicity-assessment/analyses/input/RMS_genome_cut_hg38_intervar.hg38_multianno.txt.grl_p"
input_vcf_file      <- "/Users/naqvia/d3b_codes/pathogenicity-assessment/analyses/input/RMS_genome_cut_hg38.vcf"

## filters for gnomAD
filter_gnomad_var    <- opt$gnomad_variable
filter_variant_depth <- opt$variant_depth
filter_variant_af    <- opt$variant_af

## or testing only
filter_variant_depth = 15

## output files
output_tab_file <- file.path(analysis_dir, "annotations_report.tsv") 
output_tab_abr_file  <- file.path(analysis_dir, "annotations_report.abridged.tsv")

## make vcf dataframe and add vcf_if column 
vcf_df  <-  vroom(input_vcf_file, comment = "#",delim="\t", col_names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"), show_col_types = FALSE) 
            %>% mutate(vcf_id= str_remove_all(paste (CHROM,"-",POS,"-",REF,"-",ALT), " ")) ## add vcf id column

## add clinvar table to this (INFO)
clinvar_anno_vcf_df  <-  vroom(input_clinVar_file, comment = "#", delim="\t", col_names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"), show_col_types = FALSE) %>%
                        #add c_id column
                        mutate(c_id = str_match(INFO, "MedGen:CN(\\d+);")[, 2]) %>% 
  
                         #add vcf id colum  
                         mutate(vcf_id= str_remove_all(paste ("chr",CHROM,"-",POS,"-",REF,"-",ALT), " ")) %>% 
                         full_join(vcf_df, by="vcf_id") %>% 
                         
                         #apply DP filters
                         mutate(variant_depth = if_else( as.integer( str_match(INFO.y, "DP\\=(\\d+)")[, 2])  > filter_variant_depth, "PASS","FAIL")) %>%   ## apply DP filters 
                        
                         #apply variant_af filter
                         mutate(variant_af    = if_else(as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3]) / ( (as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,2]) ) + as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3] )) > filter_variant_af, "PASS", "FAIL")) %>%
                        
                         #add star annotations to clinVar results table based on filters // ## default version
                         mutate(Stars = ifelse(grepl('CLNREVSTAT\\=criteria_provided,_single_submitter', INFO.x), "1",
                                              ifelse(grepl('CLNREVSTAT\\=criteria_provided,_multiple_submitters', INFO.x), "2",
                                                     ifelse(grepl('CLNREVSTAT\\=reviewed_by_expert_panel', INFO.x), "3",
                                                            ifelse(grepl('CLNREVSTAT\\=practice_guideline', INFO.x), "4",
                                                                   ifelse(grepl('CLNREVSTAT\\=criteria_provided,_conflicting_interpretations', INFO.x), "Needs Resolving", "0")
                                                            )))),
                               ## extract the calls and put in own column
                               final_call = str_match(INFO.x, "CLNSIG\\=(\\w+)\\;")[, 2])

## add intervar table
clinvar_anno_intervar_vcf_df  <-  vroom(input_intervar_file, comment = "#", delim="\t", show_col_types = FALSE) %>% 
                                  
                                  #add vcf id column
                                  mutate(vcf_id= str_remove_all(paste (Chr,"-",Start,"-",Ref,"-",Alt), " ")) %>% 
                                  
                                  #apply gnomad_af filters  
                                  mutate(gnomad_af = if_else( as.numeric(gnomAD_genome_ALL) > filter_variant_af, "PASS", "FAIL" )) %>%
                                  
                                  ## merge dataframe with clinvar_anno_vcf_df above
                                  full_join(clinvar_anno_vcf_df, by="vcf_id") 
  
                                  

## retrieve and store clinVar input file into table data.table::fread()
input_submissions_file_path = file.path(input_dir, "submission_summary.txt")

submission_info_df  <-  vroom(input_submissions_file_path, comment = "#",delim="\t", 
                               col_names = c("VariationID","ClinicalSignificance","DateLastEvaluated","Description","SubmittedPhenotypeInfo","ReportedPhenotypeInfo",
                                             "ReviewStatus","CollectionMethod","OriginCounts","Submitter","SCV","SubmittedGeneSymbol","ExplanationOfInterpretation"), 
                               show_col_types = FALSE) %>% mutate(c_id = str_match(ReportedPhenotypeInfo, "C(\\d+):")[, 2]) ## add c_id column to match clinvar vcf

## filter only those variants that need consensus call and find final call in submission table
entries_for_cc <- clinvar_anno_intervar_vcf_df %>% filter(Stars == "Needs Resolving", na.rm = TRUE)  
entries_for_cc <- entries_for_cc %>% left_join(submission_info_tab, by="c_id") %>% mutate(final_call=submission_info_tab$ClinicalSignificance)

##autopvs1 run
