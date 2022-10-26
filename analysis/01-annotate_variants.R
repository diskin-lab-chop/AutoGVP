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
})

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses")

input_dir   <- file.path(root_dir, "input")

## retrieve and store clinVar input file into table
input_clinVar_file <- (paste0(input_dir, "/", "clinvar_dummy.nih.norm.annot.chr1.vcf"))
clinVar_results  <-  read.delim(input_clinVar_file, comment.char = "#", sep = "\t", header=FALSE) 

##add back original column headers
colnames(clinVar_results) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

## add column "vcf_id" to clinVar results in order to cross-reference with intervar and autopvs1 table
clinvar_results <- clinVar_results %>% mutate(vcf_id= str_remove_all(paste (CHROM,"-",POS,"-",REF,"-",ALT), " "))

## add star annotations to clinVar results table based on filters // ##default version
clinvar_results <- mutate(clinvar_results, 
                                Stars = ifelse(grepl('CLNREVSTAT\\=criteria_provided,_single_submitter', INFO), "1",
                                        ifelse(grepl('CLNREVSTAT\\=criteria_provided,_multiple_submitters', INFO), "2",
                                        ifelse(grepl('CLNREVSTAT\\=reviewed_by_expert_panel', INFO), "3",
                                        ifelse(grepl('CLNREVSTAT\\=practice_guideline', INFO), "4",
                                        ifelse(grepl('CLNREVSTAT\\=criteria_provided,_conflicting_interpretations', INFO), "Needs Resolving", "0")
                                        ))))) %>% 
                          ## extract the calls and put in own column
                          mutate(Call = str_match(INFO, "CLNSIG\\=(\\w+)\\;")[, 2])

## filter only those variants that need consensus call
run_cc <- filter(clinVar_results_stars, Stars == "Needs Resolving", na.rm = TRUE)

##One Star cases that are “criteria_provided,_single_submitter” that do NOT have the B, LB, P, LP call ## need to implement this
 
## filter only those variants that need an InterVar run (No Star)
run_interVar <- filter(clinVar_results_stars, Stars == "0", na.rm = TRUE)

## retrieve and store interVar output file into table
input_intervar_file <- (paste0(input_dir, "/", "706d8be0-a723-4205-9920-a5954efc6793.hg38_multianno.txt.intervar.chr1"))
intervar_results    <-  read.delim(input_intervar_file, comment.char = "#", sep = "\t", header=FALSE) 

## autopvs1 results and file
input_autopvs1_file <- (paste0(input_dir, "/", "clinvar_dummy.nih.autopvs1.chr1.tsv"))
autopvs1_results    <-  read.delim(input_autopvs1_file, comment.char = "#", sep = "\t", header=FALSE) 

## add back headers for each table -- intervar and autopvs1 results
colnames(intervar_results) <- c("Chr","Start","End","Ref","Alt","Ref.Gene","Func.refGene",
                                "ExonicFunc.refGene","Gene.ensGene","avsnp147","AAChange.ensGene","AAChange.refGene","clinvar: Clinvar","InterVar: InterVar and Evidence", 	
                                "Freq_gnomAD_genome_ALL","Freq_esp6500siv2_all","Freq_1000g2015aug_all","CADD_raw","CADD_phred","SIFT_score","GERP++_RS","phyloP46way_placental",	
                                "dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE","Interpro_domain	AAChange.knownGene","rmsk","MetaSVM_score","Freq_gnomAD_genome_POPs",
                                "OMIM","Phenotype_MIM","OrphaNumber","Orpha	Otherinfo")

colnames(autopvs1_results) <- c("vcf_id","SYMBOL","Feature","trans_name","consequence","strength_raw","strength","criterion")

## add column "vcf_id" to intervar resultsin order to match with autopvs1 and clinvar table 
intervar_results <- intervar_results %>% mutate(vcf_id= str_remove_all(paste (Chr,"-",Start,"-",Ref,"-",Alt), " "))

## join two tables together based on variant id
combined_tab     <- inner_join(intervar_results, autopvs1_results, by="vcf_id")

## criteria to check intervar/autopvs1 to re-calculate and create a score column that will inform the new re-calculated call

#if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
mutate(combined_tab, score=ifelse( (criterion =="NF1")  |
                                   (criterion =="SS1")  |
                                   (criterion =="DEL1") |
                                   (criterion =="DEL2") |
                                   (criterion =="DUP1") |
                                   (criterion =="IC1"), "PVS1=1",""))

#if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
mutate(combined_tab, score=ifelse((criterion =="NF3")  |
                                  (criterion =="NF5")  |
                                  (criterion =="SS3") |
                                  (criterion =="SS5") |
                                  (criterion =="SS8") |
                                  (criterion =="SS10")|
                                  (criterion =="DEL8") |
                                  (criterion =="DEL6") |
                                  (criterion =="DEL10")|
                                  (criterion =="DUP3")|
                                  (criterion =="IC2"), paste("PVS1=0","PS=PS+1"), ""))

#if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
mutate(combined_tab, score=ifelse((criterion =="NF6")  |
                                  (criterion =="SS6")  |
                                  (criterion =="SS9") |
                                  (criterion =="DEL7") |
                                  (criterion =="DEL11") |
                                  (criterion =="IC3"), paste("PVS1=1","PM=PM+1"),""))

#if criterion is IC4 then PVS1 = 0; PP = PP+1;
mutate(combined_tab, score=ifelse((criterion =="IC4"), paste("PVS1=0","PP=PP+1"),""))

#if criterion is na then PVS1 = 0;
mutate(combined_tab, score=ifelse((criterion =="na"), "PVS1=0",""))

## add new call based on new scoring metric 




