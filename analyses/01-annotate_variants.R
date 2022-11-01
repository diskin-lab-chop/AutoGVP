################################################################################
# 01-annotate_variants.R
# written by Ammar Naqvi
#
# This script annotates variants based on clinVar, and recomputes score and call 
# for intervar and autoPVS1
#
# usage: Rscript 01-annotate_variants.R--vcf <vcf file> --intervar <intervar file> --autopvs1 <autopvs1 file>
#
################################################################################

suppressPackageStartupMessages({
  library("tidyverse")
  library("optparse")
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
  make_option(c("--autopvs1"), type = "character",
              help = "input autopvs1 file")
)

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters
input_clinVar_file <-  opt$vcf
input_intervar_file <- opt$intervar
input_autopvs1_file <- opt$autopvs1

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
  
  ## add column for individual scores that will be re-calculated if we need to adjust using autoPVS1 result
  mutate(evidencePVS1 = str_match(`InterVar: InterVar and Evidence`, "PVS1\\=(\\d+)\\s")[, 2], 
         evidencePS = str_match(`InterVar: InterVar and Evidence`, "PS\\=(\\d+)\\s")[, 2], 
         evidencePM = str_match(`InterVar: InterVar and Evidence`, "PM\\=(\\d+)\\s")[, 2], 
         evidencePP = str_match(`InterVar: InterVar and Evidence`, "PP\\=(\\d+)\\s")[, 2]) %>% 
         replace(is.na(.), "0")

  combined_tab_for_intervar <- combined_tab_for_intervar %>% 
  
  ## if PVS1=0, take intervar call
  mutate(final_call = if_else(evidencePVS1 == 0, str_match(`InterVar: InterVar and Evidence`, "InterVar\\:\\s+(.+?(?=\\sPVS))")[, 2], "recalculate")) %>%
  
  ## criteria to check intervar/autopvs1 to re-calculate and create a score column that will inform the new re-calculated final call 
  #if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
  mutate(evidencePVS1 = if_else( (criterion == "NF1" | criterion == "SS1" |
                                   criterion == "DEL1" | criterion == "DEL2" |
                                   criterion == "DUP1" | criterion == "IC1") & evidencePVS1 == 1, "1", evidencePVS1)) %>% 

  #if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
  mutate(evidencePVS1 = if_else(  (criterion == "NF3"  | criterion == "NF5" |
                                   criterion == "SS3" | criterion == "SS5" |
                                   criterion == "SS8" | criterion =="SS10" |
                                   criterion =="DEL8" | criterion =="DEL6" |
                                   criterion =="DEL10" | criterion =="DUP3" |
                                   criterion =="IC2") & evidencePVS1 == 1, "0", evidencePVS1)) %>%
  mutate(evidencePS =  if_else(   (criterion == "NF3"  | criterion == "NF5" |
                                  criterion == "SS3" | criterion == "SS5" |
                                  criterion == "SS8" | criterion =="SS10" |
                                  criterion =="DEL8" | criterion =="DEL6" |
                                  criterion =="DEL10" | criterion =="DUP3" |
                                  criterion =="IC2") & evidencePVS1 == 1, as.numeric(evidencePS)+1, as.double(evidencePS))) %>%
          
  #if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
  mutate(evidencePVS1 = if_else( (criterion == "NF6" | criterion == "SS6" |
                                  criterion == "SS9" | criterion == "DEL7" |
                                  criterion == "DEL11" | criterion == "IC3") & evidencePVS1 == 1, "0", evidencePVS1)) %>% 
  mutate(evidencePM =   if_else( (criterion == "NF6" | criterion == "SS6" |
                                  criterion == "SS9" | criterion == "DEL7" |
                                  criterion == "DEL11" | criterion == "IC3") & evidencePVS1 == 1, as.numeric(evidencePM)+1, as.double(evidencePM))) %>% 
  
  #if criterion is IC4 then PVS1 = 0; PP = PP+1;
  mutate(evidencePVS1 = if_else((criterion == "IC4") & evidencePVS1 == 1, "0", evidencePVS1)) %>%
  mutate(evidencePP   = if_else((criterion == "IC4") & evidencePVS1 == 1, as.numeric(evidencePP)+1, as.double(evidencePP)) %>%
    
  #if criterion is na then PVS1 = 0;
  mutate(evidencePVS1 = if_else( (criterion == "na") & evidencePVS1 == 1, 0, evidencePVS1))

## add new call based on new re-calculated scores 