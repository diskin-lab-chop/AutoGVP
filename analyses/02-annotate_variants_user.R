################################################################################
# 02-annotate_variants_user.R
# written by Ammar Naqvi
#
# This script annotates variants based on clinVar, and recomputes score and calls 
# for intervar and autoPVS1
#
# usage: Rscript 02-annotate_variants_user.R --vcf <VEP annotated vcf file> 
#                                       --intervar <intervar file> 
#                                       --autopvs1 <autopvs1 file>
#                                       --clinvar  <e.g. clinvar_20211225.vcf.gz>
#                                       --submission <variant_summary file>
#                                       --submission_summary <submission_summary file>
#                                       --gnomad_variable 'gnomad_3_1_1_AF_non_cancer'
#                                       --gnomad_af <numeric>
#                                       --variant_depth <integer>
#                                       --variant_af <numeric>
#                                       --summary_level_vcf <numeric>
#                                       --output <string>
################################################################################

## load libraries 
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
  make_option(c("--multianno"), type = "character",
              help = "input multianno file"),
  make_option(c("--autopvs1"), type = "character",
              help = "input autopvs1 file"),
  make_option(c("--clinvar"), type = "character",
              help = "specific clinVar file (format: clinvar_20211225.vcf.gz)"), ## https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022clinvar_20211225.vcf.gz 
  make_option(c("--submission"), type = "character",
              help = "specific submission file (format: variant_summary_2023-02.txt.gz)"),
  make_option(c("--submission_summary"), type = "character",
              help = "specific submission summary file (format: submission_summary.txt)"),
  make_option(c("--gnomad_variable"), type = "character",default = "Freq_gnomAD_genome_ALL",
              help = "gnomAD variable"),
  make_option(c("--gnomad_af"), type = "numeric", default = 0.001,
              help = "genomAD AF filter"),
  make_option(c("--variant_depth"), type = "integer",default = 15,
              help = "variant depth filter"),
  make_option(c("--variant_af"), type = "numeric", default = .2,
              help = "variant AF cut-off"),
  make_option(c("--summary_level"), type = "character", default = "F",
              help = "summary_level T/F"),
  make_option(c("--output"), type = "character", default = "out",
              help = "output name")) 


opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (read)
input_clinVar_file  <-  opt$clinvar
input_intervar_file <- opt$intervar
input_autopvs1_file <- opt$autopvs1
input_vcf_file <- opt$vcf
sample_name <- opt$sample_name
input_submission_file  <-  opt$submission
input_summary_submission_file <- opt$submission_summary
input_multianno_file <- opt$multianno
summary_level <- opt$summary_level
output_name   <- opt$output

## filters for gnomAD
filter_gnomad_var    <- opt$gnomad_variable
filter_variant_depth <- opt$variant_depth
filter_variant_af    <- opt$variant_af

## for testing purposes only ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
input_clinVar_file  <- file.path(input_dir,"clinvar.vcf.gz")
input_intervar_file <- file.path(input_dir,"RMS_hg38_selected_InterVar.hg38_multianno.txt.intervar")
input_multianno_file <- file.path(input_dir, "RMS_selected.hg38_multianno.txt")
input_vcf_file      <- file.path(input_dir,"RMS_hg38_selected_VEP_anntoated.vcf")
input_autopvs1_file <- file.path(input_dir,"RMS_hg38_selected_autopvs1.txt")
input_submission_file <- file.path(input_dir,"variant_summary.txt")
input_summary_submission_file <- file.path(input_dir,"submission_summary.txt")

sample_name <-  "test"
gnomad_variable <- "Freq_gnomAD_genome_ALL"
gnomad_af <- 0.001
filter_variant_depth = 15
## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

## output files
output_tab_file      <- file.path(analysis_dir, paste0(output_name,".annotations_report.tsv")) 
output_tab_abr_file  <- file.path(analysis_dir, paste0(output_name,".annotations_report.abridged.tsv"))
output_tab_dev_file  <- file.path(analysis_dir, paste0(output_name,".annotations_report.abridged.dev.tsv"))

## allocate more memory capacity
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

## make vcf dataframe and add vcf_if column 
vcf_df  <-  vroom(input_vcf_file, comment = "#",delim="\t", col_names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"), trim_ws = TRUE, show_col_types = FALSE) %>%
            mutate(vcf_id= str_remove_all(paste (CHROM,"-",POS,"-",REF,"-",ALT), " ")) %>% 
            mutate(vcf_id = str_replace_all(vcf_id, "chr", ""))  ## add vcf id column

## add clinvar table to this (INFO)
clinvar_anno_vcf_df  <- vroom(input_clinVar_file, comment = "#", delim="\t", col_names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"),trim_ws = TRUE, show_col_types = FALSE) %>%

                         #add vcf id column  
                         mutate(vcf_id= str_remove_all(paste (CHROM,"-",POS,"-",REF,"-",ALT), " ")) %>% 
                         semi_join(vcf_df, by="vcf_id") %>%  
                         mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) %>% 

                         #add star annotations to clinVar results table based on filters // ## default version
                         mutate(Stars = ifelse(grepl('CLNREVSTAT\\=criteria_provided,_single_submitter', INFO), "1",
                                              ifelse(grepl('CLNREVSTAT\\=criteria_provided,_multiple_submitters', INFO), "2",
                                                     ifelse(grepl('CLNREVSTAT\\=reviewed_by_expert_panel', INFO), "3",
                                                            ifelse(grepl('CLNREVSTAT\\=practice_guideline', INFO), "4",
                                                                   ifelse(grepl('CLNREVSTAT\\=criteria_provided,_conflicting_interpretations', INFO), "1NR", "0")
                                                            )))),
                               ## extract the calls and put in own column
                               #final_call = str_match(INFO, "CLNSIG\\=(\\w+([\\/\\|]\\w+)*)\\;")[, 2])
                               final_call = str_match(INFO, "CLNSIG\\=(\\w+)(\\|]\\w+)*\\;")[, 2])


## if conflicting intrep. take the call with most calls in CLNSIGCONF field
for(i in 1:nrow(clinvar_anno_vcf_df)) {
  entry <- clinvar_anno_vcf_df[i,]
  if(entry$Stars != "1NR")
  {
    next;
  }
  
  conf_section <-str_match(entry$INFO, "CLNSIGCONF\\=.+\\;CLNVC")  ## part to parse and count calls
  call_names <- c("Pathogenic","Likely_pathogenic","Benign","Likely_benign","Uncertain_significance")

  P  <-  (str_match(conf_section, "Pathogenic\\((\\d+)\\)")[,2])
  LP <-  (str_match(conf_section, "Likely_pathogenic\\((\\d+)\\)")[,2])
  B  <-  (str_match(conf_section, "Benign\\((\\d+)\\)")[,2])
  LB <-  (str_match(conf_section, "Likely_benign\\((\\d+)\\)")[,2])
  U  <-  (str_match(conf_section, "Uncertain_significance\\((\\d+)\\)")[,2])

  ## make vector out of possible calls to get max 
  calls          <- c(P,LP,B,LB,U)

  if ( length( which( calls == max(calls,na.rm = TRUE) ) ) > 1 )
  {
    next;
  }

  highest_ind    <- which.max(calls)
  consensus_call <- call_names[highest_ind]
  
  clinvar_anno_vcf_df[i,]$final_call = consensus_call
}

## filters (gnomAD, variant_depth, variant AF if applicable)
if(summary_level == "T")
{
    clinvar_anno_vcf_df %>%  if_else( as.integer( str_match(INFO, "DP\\=(\\d+)")[, 2])  > filter_variant_depth, "PASS","FAIL") %>% 
    mutate(gnomad_af     = if_else( as.numeric( str_match(INFO, "gnomad_3_1_1_AF_non_cancer\\=(0\\.\\d+)\\;")[,2])  > filter_variant_af, "PASS","FAIL")) %>% 
    mutate(variant_af    = if_else(as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3]) / ( (as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,2]) ) + as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3] )) > filter_variant_af, "PASS", "FAIL"))
} 

## store variants without clinvar info
clinvar_anti_join_vcf_df  <- anti_join(vcf_df, clinvar_anno_vcf_df, by="vcf_id") %>%
                             mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) %>% 
                             mutate(CHROM = str_replace_all(CHROM, "chr", ""))   %>% 
                             dplyr::rename(rs_id = ID)

## get latest calls from submission files 
submission_summary_df <- vroom(input_summary_submission_file, comment = "#",delim="\t", 
                            col_names = c("VariationID","ClinicalSignificance","DateLastEvaluated",       
                                          "Description","SubmittedPhenotypeInfo","ReportedPhenotypeInfo",   
                                          "ReviewStatus","CollectionMethod","OriginCounts","Submitter",
                                          "SCV","SubmittedGeneSymbol","ExplanationOfInterpretation"),
                            show_col_types = FALSE) %>% dplyr::select("VariationID","ClinicalSignificance") %>% 
                            group_by(VariationID) %>%
                            arrange(ClinicalSignificance) %>%
                            dplyr::slice(1) %>%
                            ungroup


submission_info_df  <-  vroom(input_submission_file, delim="\t",
                              col_types = c(ReferenceAlleleVCF = "c",AlternateAlleleVCF= "c",PositionVCF="i",VariationID="n" ),
                                 show_col_types = FALSE) %>% 
  
                             #add vcf id column  
                             mutate(vcf_id= str_remove_all(paste (Chromosome,"-",PositionVCF,"-",ReferenceAlleleVCF,"-",AlternateAlleleVCF), " ")) %>% 
                             mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) %>% 
                             mutate(VariationID=as.double(noquote(VariationID)))
                                
## join submission files to ensure we have vcf id to match with other tables
submission_final_df <- inner_join(submission_summary_df,submission_info_df, by="VariationID" )

## filter only those variants that need consensus call and find  call in submission table
entries_for_cc <-  filter(clinvar_anno_vcf_df, Stars == "1NR", final_call !="Benign",final_call !="Pathogenic", final_call != "Likely_benign",final_call !="Likely_pathogenic", final_call != "Uncertain_significance")
entries_for_cc_in_submission <- inner_join(submission_final_df,entries_for_cc, by="vcf_id") %>% dplyr::mutate(final_call=ClinicalSignificance.x) %>% 
                                dplyr::select(vcf_id,ClinicalSignificance.x,final_call) %>% dplyr::rename("ClinicalSignificance"=ClinicalSignificance.x)
                 
## one Star cases that are “criteria_provided,_single_submitter” that do NOT have the B, LB, P, LP, VUS call must also go to intervar
#additional_intervar_cases <-  filter(clinvar_anno_vcf_df, Stars == "1", final_call!="Benign",final_call!="Pathogenic", final_call != "Likely_benign",final_call!="Likely_pathogenic", final_call != "Uncertain_significance")
additional_intervar_cases <-  filter(clinvar_anno_vcf_df, final_call!="Benign",final_call!="Pathogenic", final_call != "Likely_benign",final_call!="Likely_pathogenic", final_call != "Uncertain_significance")


clinvar_anti_join_vcf_df <- clinvar_anti_join_vcf_df %>% mutate(QUAL = as.character(QUAL))

## filter only those variant entries that need an InterVar run (No Star) and add the additional intervar cases from above
entries_for_intervar <- filter(clinvar_anno_vcf_df, Stars == "0", na.rm = TRUE) %>% bind_rows((additional_intervar_cases)) %>% bind_rows(clinvar_anti_join_vcf_df)

## get vcf ids that need intervar run
vcf_to_run_intervar <- entries_for_intervar$vcf_id

## get multianno file to add by correct vcf_id
multianno_df  <-  vroom(input_multianno_file, delim="\t",trim_ws = TRUE, col_names = TRUE, show_col_types = FALSE) %>% 
                    mutate(vcf_id= str_remove_all(paste (Chr,"-",Otherinfo5,"-",Otherinfo7,"-",Otherinfo8), " ")) %>% 
                    mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) %>% 
                    group_by(vcf_id) %>%
                    arrange(vcf_id) %>%
                    filter(row_number()==1) %>% 
                    ungroup

## add intervar table
clinvar_anno_intervar_vcf_df  <-  vroom(input_intervar_file, delim="\t",trim_ws = TRUE, col_names = TRUE, show_col_types = TRUE) %>% 
                                  #slice(-1) %>% 
                                  mutate(var_id= str_remove_all(paste (`#Chr`,"-",Start,"-",End,"-",Ref,"-",Alt), " ")) %>% 
                                  group_by(var_id) %>%
                                  arrange(var_id) %>%
                                  filter(row_number()==1) %>% 
                                  ungroup

## exit if the total number of variants differ in these two tables to ensure we annotate with the correct vcf so we can match back to clinVar and other tables
if( tally(multianno_df) != tally(clinvar_anno_intervar_vcf_df) ) {
  stop("intervar and multianno files of diff lengths") 
}

## combine the intervar and multianno tables by the appropriate vcf id
clinvar_anno_intervar_vcf_df <- mutate(multianno_df,clinvar_anno_intervar_vcf_df )  %>% dplyr::select(vcf_id,`InterVar: InterVar and Evidence`, Ref.Gene,Func.refGene,ExonicFunc.refGene)

## populate consensus call variants with invervar info
entries_for_cc_in_submission_w_intervar <- inner_join(clinvar_anno_intervar_vcf_df,entries_for_cc_in_submission, by="vcf_id") %>% 
                                            dplyr::select(vcf_id,`InterVar: InterVar and Evidence`,Ref.Gene,Func.refGene,ExonicFunc.refGene) %>% 
                                            dplyr::rename("Intervar_evidence"=`InterVar: InterVar and Evidence`)


clinvar_anno_intervar_vcf_df <- clinvar_anno_intervar_vcf_df %>%  anti_join(entries_for_cc_in_submission, by="vcf_id") %>%
  ## add column for individual scores that will be re-calculated if we need to adjust using autoPVS1 result
  ## note: ignore PP5 score
  mutate(evidencePVS1 = str_match(`InterVar: InterVar and Evidence`, "PVS1\\=(\\d+)\\s")[, 2]) %>%
  mutate(evidenceBA1 = str_match(`InterVar: InterVar and Evidence`, "BA1\\=(\\d+)\\s")[, 2]) %>%
  mutate( evidencePS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidencePM = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPM\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidencePP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))[-5])) ) %>%
  mutate( evidenceBS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidenceBP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  
  ## merge dataframe with clinvar_anno_vcf_df above
  full_join(clinvar_anno_vcf_df, by="vcf_id") 

## autopvs1 results
autopvs1_results    <-  read_tsv(input_autopvs1_file, col_names = TRUE) %>%
  mutate(vcf_id = str_remove_all(paste (vcf_id), " ")) %>% 
  mutate(vcf_id = str_replace_all(vcf_id, "chr", "")) 

## join all three tables together based on variant id that need intervar run
combined_tab_for_intervar <- autopvs1_results %>%
  inner_join(clinvar_anno_intervar_vcf_df, by="vcf_id") %>% 
  dplyr::filter(vcf_id %in% entries_for_intervar$vcf_id) %>% dplyr::select(vcf_id,`InterVar: InterVar and Evidence`, criterion, evidencePVS1, evidenceBA1, evidencePS, evidencePM, evidencePP, evidenceBS, evidenceBP) %>% 

  ## indicate if recalculated 
  mutate(intervar_adjusted_call = if_else( (evidencePVS1 == 0), "No", "Yes")) %>%  
  
  ## criteria to check intervar/autopvs1 to re-calculate and create a score column that will inform the new re-calculated final call
  #if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
  mutate(evidencePVS1 = if_else( (criterion == "NF1" | criterion == "SS1" |
                                    criterion == "DEL1" | criterion == "DEL2" |
                                    criterion == "DUP1" | criterion == "IC1") & evidencePVS1 == 1, "1", evidencePVS1)) %>%
  
  #if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
  mutate(evidencePS =  if_else(   (criterion == "NF3"  | criterion == "NF5" |
                                     criterion == "SS3" | criterion == "SS5" |
                                     criterion == "SS8" | criterion =="SS10" |
                                     criterion =="DEL8" | criterion =="DEL6" |
                                     criterion =="DEL10" | criterion =="DUP3" |
                                     criterion =="IC2") & evidencePVS1 == 1, as.numeric(evidencePS)+1, as.double(evidencePS))) %>%
  
  mutate(evidencePVS1 = if_else(  (criterion == "NF3"  | criterion == "NF5" |
                                     criterion == "SS3" | criterion == "SS5" |
                                     criterion == "SS8" | criterion =="SS10" |
                                     criterion =="DEL8" | criterion =="DEL6" |
                                     criterion =="DEL10" | criterion =="DUP3" |
                                     criterion =="IC2") & evidencePVS1 == 1, "0", evidencePVS1)) %>%
  
  #if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
  mutate(evidencePM =   if_else( (criterion == "NF6" | criterion == "SS6" |
                                    criterion == "SS9" | criterion == "DEL7" |
                                    criterion == "DEL11" | criterion == "IC3") & evidencePVS1 == 1, as.numeric(evidencePM)+1, as.double(evidencePM))) %>%
  
  mutate(evidencePVS1 = if_else( (criterion == "NF6" | criterion == "SS6" |
                                    criterion == "SS9" | criterion == "DEL7" |
                                    criterion == "DEL11" | criterion == "IC3") & evidencePVS1 == 1, "0", evidencePVS1)) %>%
  
  #if criterion is IC4 then PVS1 = 0; PP = PP+1;
  mutate(evidencePP   = if_else((criterion == "IC4") & evidencePVS1 == 1, as.numeric(evidencePP)+1, as.double(evidencePP))) %>%
  mutate(evidencePVS1 = if_else((criterion == "IC4") & evidencePVS1 == 1, "0", evidencePVS1)) %>%
  
  #if criterion is na then PVS1 = 0;
  mutate(evidencePVS1 = if_else( (criterion == "na") & evidencePVS1 == 1, 0, as.double(evidencePVS1))) %>%
  
  ## adjust variables based on given rules described in README
  mutate(final_call = ifelse( (evidencePVS1)   == 1 &
                              ((evidencePS   >= 1) |
                               (evidencePM   >=2 ) |
                               (evidencePM   >= 1 & evidencePP ==1) |
                               (evidencePP   >=2 ) ), "Pathogenic",
                            ifelse( (evidencePVS1   == 1 & evidencePS >= 2), "Pathogenic",
                                    ifelse( (evidencePVS1)   == 1 & (evidencePS == 1 &
                                            (evidencePM   >= 3) |
                                            (evidencePM   ==2 & evidencePP >=2 ) |
                                            (evidencePM == 1 & evidencePP >=4 ) ), "Pathogenic",
                                            ifelse( (evidencePVS1 == 1 & evidencePM == 1) |
                                                    (evidencePS==1 & evidencePM >= 1) |
                                                    (evidencePS==1 & evidencePP >=2 ) |
                                                    (evidencePM >= 3) |
                                                    (evidencePM ==2 & evidencePP>=2 ) |
                                                    (evidencePM == 1 & evidencePP>=4), "Likely_pathogenic",
                                                    ifelse( (evidenceBA1 == 1) |
                                                            (evidenceBS   >= 2), "Benign",
                                                            ifelse( (evidenceBS == 1 & evidenceBP == 1) |
                                                                    (evidenceBP   >= 2), "Likely_benign",  
                                                                    ifelse( evidencePVS1 == 0, str_match(`InterVar: InterVar and Evidence`, "InterVar\\:\\s+(.+?(?=\\sPVS))")[, 2],"Uncertain_significance"))))))))





## merge tables together (clinvar and intervar) and write to file
master_tab <- full_join(clinvar_anno_intervar_vcf_df,combined_tab_for_intervar, by="vcf_id" ) 
master_tab <- master_tab %>% mutate(intervar_adjusted_call = coalesce(intervar_adjusted_call, "Not adjusted, clinVar")) %>% 
                             mutate(evidencePVS1 = coalesce(as.double(evidencePVS1.x, evidencePVS1.y) )) %>%
                             mutate(evidenceBA1 = coalesce(as.double(evidenceBA1.x, evidenceBA1.y) )) %>% 
                             mutate(evidencePS = coalesce(as.double(evidencePS.x, evidencePS.y) )) %>% 
                             mutate(evidencePM = coalesce(as.double(evidencePM.x, evidencePM.y) )) %>% 
                             mutate(evidencePP = coalesce(as.double(evidencePP.x, evidencePP.y) )) %>% 
                             mutate(evidenceBS = coalesce(as.double(evidenceBS.x, evidenceBS.y) )) %>% 
                             mutate(evidenceBP = coalesce(as.double(evidenceBP.x, evidenceBP.y) )) %>% 
                             mutate(Intervar_evidence = coalesce(`InterVar: InterVar and Evidence.x`, `InterVar: InterVar and Evidence.y`)) %>% 
  
                             # replace second final call with the second one because we did not use interVar results
                             mutate(final_call.x = if_else(intervar_adjusted_call=="No" & Stars=="0", final_call.y, final_call.x))  

## combine final calls into one choosing the appropriate final call                             
master_tab <- master_tab %>% mutate(final_call = coalesce(final_call.x, final_call.y))

## remove older columns
master_tab <- master_tab %>% dplyr::select(-c (evidencePVS1.x,evidencePVS1.y,evidenceBA1.x, evidenceBA1.y, evidencePS.x, evidencePS.y, evidencePM.x, evidencePM.y, evidencePP.x, evidencePP.y, evidenceBS.x, evidenceBS.y, evidenceBP.x, evidenceBP.y,
                                    `InterVar: InterVar and Evidence.x`, `InterVar: InterVar and Evidence.y`,final_call.x,final_call.y))

## reformat columns
master_tab  <- full_join(master_tab,entries_for_cc_in_submission, by="vcf_id") %>% 
               mutate(final_call = coalesce(final_call.y, final_call.x)) %>% 
               full_join(entries_for_cc_in_submission_w_intervar, by="vcf_id") %>%
               mutate(Intervar_evidence = coalesce(Intervar_evidence.y, Intervar_evidence.x)) %>%
               mutate(Ref.Gene = coalesce(Ref.Gene.y, Ref.Gene.x)) 
  


## write to output files
# master version with all columns
# write.table(master_tab, output_tab_file, append = FALSE, sep = "\t", dec = ".",row.names = FALSE, quote = FALSE, col.names = TRUE)

# development version 
#results_tab_dev <- master_tab %>% dplyr::select(vcf_id, Stars, Intervar_evidence, evidencePVS1,evidencePS,evidencePM,evidencePP,evidencePP, evidenceBA1, evidenceBS, evidenceBP, intervar_adjusted_call, final_call)
#write.table(results_tab_dev, output_tab_dev_file, append = FALSE, sep = "\t", dec = ".",
#            row.names = FALSE, quote = FALSE, col.names = TRUE)

# abridged version
results_tab_abridged <- master_tab %>% dplyr::select(vcf_id, Ref.Gene, Stars, Intervar_evidence,intervar_adjusted_call,final_call)

## if conflicting intrep. take the call with most calls in CLNSIGCONF field
for(i in 1:nrow(results_tab_abridged)) {
  entry <- results_tab_abridged[i,]
  if(is.na(entry$final_call) || (entry$final_call !="Pathogenic" && 
     entry$final_call != "Likely_benign" &&  entry$final_call !="Likely_pathogenic"
     && entry$final_call != "Uncertain_significance"  &&  entry$final_call !="Benign"  
     &&  entry$final_call !="Uncertain significance"  &&  entry$final_call !="Likely benign") )
  {
    print (results_tab_abridged[i,]$final_call)
    
    new_call <- str_match(results_tab_abridged[i,]$Intervar_evidence, "InterVar\\:\\s(\\w+\\s\\w+)*")[,2]
    results_tab_abridged[i,]$final_call = new_call
  }
}

results_tab_abridged <- results_tab_abridged %>% mutate(final_call = replace(final_call, final_call == "Likely benign", "Likely_benign"))
results_tab_abridged <- results_tab_abridged %>% mutate(final_call = replace(final_call, final_call == "Uncertain significance", "Uncertain_significance"))
results_tab_abridged <- results_tab_abridged %>% mutate(final_call = replace(final_call, final_call == "Benign PVS1", "Benign"))
results_tab_abridged <- results_tab_abridged %>% mutate(final_call = replace(final_call, final_call == "Pathogenic PVS1","Pathogenic"))
results_tab_abridged <- results_tab_abridged %>% mutate(final_call = replace(final_call, final_call == "Likely pathogenic","Likely_pathogenic"))




write.table(results_tab_abridged, output_tab_abr_file, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, quote = FALSE, col.names = TRUE)
