################################################################################
# 01-annotate_variants.R
# written by Ammar Naqvi
#
# This script annotates variants based on clinVar, and recomputes score and call 
# for intervar and autoPVS1
#
# usage: Rscript 01-annotate_variants.R
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
  
})

# Get `magrittr` pipe
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
  make_option(c("--autopvs1"), type = "character",
              help = "input autopvs1 file")
)

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters
input_clinVar_file <-  opt$vcf
input_intervar_file <- opt$intervar
input_autopvs1_file <- opt$autopvs1

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses")

## retrieve and store clinVar input file into table
clinVar_results  <-  read_tsv(input_clinVar_file, comment = "#",
                              col_names = c("CHROM","START","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"))

## add column "vcf_id" to clinVar results in order to cross-reference with intervar and autopvs1 table
clinvar_results <- clinVar_results %>%
  mutate(vcf_id= str_remove_all(paste (CHROM,"-",START,"-",REF,"-",ALT), " "),
         ## add star annotations to clinVar results table based on filters // ## default version
         Stars = ifelse(grepl('CLNREVSTAT\\=criteria_provided,_single_submitter', INFO), "1",
                 ifelse(grepl('CLNREVSTAT\\=criteria_provided,_multiple_submitters', INFO), "2",
                 ifelse(grepl('CLNREVSTAT\\=reviewed_by_expert_panel', INFO), "3",
                 ifelse(grepl('CLNREVSTAT\\=practice_guideline', INFO), "4",
                 ifelse(grepl('CLNREVSTAT\\=criteria_provided,_conflicting_interpretations', INFO), "Needs Resolving", "0")
                                      )))),
        
          ## extract the calls and put in own column
         Call = str_match(INFO, "CLNSIG\\=(\\w+)\\;")[, 2])

## filter only those variants that need consensus call
entries_for_cc <- clinvar_results %>%
          filter(Stars == "Needs Resolving", na.rm = TRUE)

##one Star cases that are “criteria_provided,_single_submitter” that do NOT have the B, LB, P, LP call must also go to intervar
additional_intervar_cases <- clinvar_results %>% filter(Stars == "0", Call!="Benign",Call!="Pathogenic", Call != "Likely_benign",Call!="Likeley_pathogenic")

## filter only those variant entries that need an InterVar run (No Star)
entries_for_intervar <- clinvar_results %>%
  filter(Stars == "0", na.rm = TRUE) %>% bind_rows((additional_intervar_cases))

clinvar_results %>%
  filter(Stars == "0" & Call != "Benign")

## get vcf ids that need intervar run
vcf_to_run_intervar <- entries_for_intervar$vcf_id

## retrieve and store interVar output file into table
intervar_results    <-  read_tsv(input_intervar_file, comment = "#", col_names = c("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene",
                                                                                   "ExonicFunc.refGene","Gene.ensGene","avsnp147","AAChange.ensGene","AAChange.refGene","clinvar: Clinvar","InterVar: InterVar and Evidence",
                                                                                   "Freq_gnomAD_genome_ALL","Freq_esp6500siv2_all","Freq_1000g2015aug_all","CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phyloP46way_placental",
                                                                                   "dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE","Interpro_domain	AAChange.knownGene","rmsk","MetaSVM_score","Freq_gnomAD_genome_POPs",
                                                                                   "OMIM","Phenotype_MIM","OrphaNumber","Orpha	Otherinfo"))

## autopvs1 results and file
autopvs1_results    <-  read_tsv(input_autopvs1_file, comment = "#", col_names = c("vcf_id","SYMBOL","Feature","trans_name","consequence","strength_raw","strength","criterion"))

## add column "vcf_id" to intervar resultsin order to match with autopvs1 and clinvar table
intervar_results <- intervar_results %>%
  mutate(vcf_id = str_remove_all(paste (Chr,"-",Start,"-",Ref,"-",Alt), " "))

## join all three tables together based on variant id that need intervar run
combined_tab_for_intervar <- autopvs1_results %>%
  inner_join(intervar_results, by="vcf_id") %>%
  filter(vcf_id %in% entries_for_intervar$vcf_id) %>% 

  ## criteria to check intervar/autopvs1 to re-calculate and create a score column that will inform the new re-calculated call
  #if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
  mutate(score = case_when(criterion == "NF1" | criterion == "SS1" |
                             criterion == "DEL1" | criterion == "DEL2" |
                             criterion == "DUP1" | criterion == "IC1" ~  "PVS1=1",
                           TRUE ~ NA_character_)) %>%

#if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
mutate(score = case_when(criterion == "NF3"  | criterion == "NF5" |
                           criterion == "SS3" | criterion == "SS5" |
                           criterion == "SS8" | criterion =="SS10" |
                           criterion =="DEL8" | criterion =="DEL6" |
                           criterion =="DEL10" | criterion =="DUP3" |
                           criterion =="IC2" ~ "PVS1=0 PS=PS+1",
                         TRUE ~ as.character(score))) %>%

#if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
mutate(score = case_when(criterion == "NF6" | criterion == "SS6" |
                           criterion == "SS9" | criterion == "DEL7" |
                           criterion == "DEL11" | criterion == "IC3" ~ "PVS1=1 PM=PM+1",
                         TRUE ~ as.character(score))) %>%

#if criterion is IC4 then PVS1 = 0; PP = PP+1;
mutate(score = case_when(criterion == "IC4" ~ "PVS1=0 PP=PP+1",
                         TRUE ~ as.character(score))) %>%

#if criterion is na then PVS1 = 0;
mutate(score = case_when(criterion == "na" ~ "PVS1=0",
                         TRUE ~ as.character(score)))

## add new call based on new scoring metric