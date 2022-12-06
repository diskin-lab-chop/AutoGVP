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
  make_option(c("--autopvs1"), type = "character",
              help = "input autopvs1 file"),
  make_option(c("--clinvar"), type = "character",
              help = "specific clinVar version (format: clinvar_20211225)"), ## https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022clinvar_20211225.vcf.gz 
  make_option(c("--gnomad_variable"), type = "character",default = "gnomad_3_1_1_AF_non_cancer",
              help = "gnomAD variable"),
  make_option(c("--gnomad_af"), type = "numeric", default = 0.001,
              help = "genomAD AF filter"),
  make_option(c("--variant_depth"), type = "integer",default = 15,
              help = "variant depth filter"),
  make_option(c("--variant_af"), type = "numeric", default = .2,
             help = "variant AF cut-off")) 

opt <- parse_args(OptionParser(option_list = option_list))

## get input files from parameters (reqd)
input_clinVar_file  <-  opt$vcf
input_intervar_file <- opt$intervar
input_autopvs1_file <- opt$autopvs1
clinvar_ver = opt$clinvar

 input_clinVar_file  <- "/Users/naqvia/d3b_codes/pathogenicity-assessment/analyses/input/clinvar_dummy.nih.norm.annot.chr17.vcf" 

 input_intervar_file <- "/Users/naqvia/d3b_codes/pathogenicity-assessment/analyses/input/clinvar_input_wf_test.hg38_multianno.chr17.txt.intervar"

 input_autopvs1_file <- "/Users/naqvia/d3b_codes/pathogenicity-assessment/analyses/input/clinvar_dummy.nih.autopvs1.chr17.tsv"

 ## filters for gnomAD
filter_gnomad_var    <- opt$gnomad_variable
filter_variant_depth <- opt$variant_depth
filter_variant_af    <- opt$variant_af

output_tab_file <- file.path(analysis_dir, "annotations_report.tsv") 
output_tab_abr_file  <- file.path(analysis_dir, "annotations_report.abridged.tsv")


## retrieve clinvar vcf if specified and download to input folder
if(is.character(clinvar_ver)){
  ## generate full path to download
  clinvar_ftp_path <- paste("ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/",clinvar_ver,".vcf.gz", sep="")
  system(paste("wget -l 3", clinvar_ftp_path, "-P ",input_dir), wait = TRUE)
}  

# ## function for gnomAD, variant af and depth filtering 
# gnomad_filtering <- function(clinvar_vcf) {
#   clinvar_vcf <- clinvar_vcf %>% 
#   mutate(variant_depth = if_else( as.integer( str_match(INFO, "DP\\=(\\d+)")[, 2])  > filter_variant_depth, "PASS","FAIL")) %>% 
#   mutate(gnomad_af     = if_else( as.numeric( str_match(INFO, "gnomad_3_1_1_AF_non_cancer\\=(0\\.\\d+)\\;")[,2])  > filter_variant_af, "PASS","FAIL")) %>% 
#   mutate(variant_af    = if_else(as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3]) / ( (as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,2]) ) + as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3] )) > filter_variant_af, "PASS", "FAIL"))
#   return(clinvar_vcf)
# }

## retrieve and store clinVar input file into table data.table::fread()
clinVar_results  <-  vroom(input_clinVar_file, comment = "#",delim="\t", col_names = c("CHROM","START","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"), show_col_types = FALSE)

## filter for gnomad, read depth and AF
#clinVar_results <- gnomad_filtering(clinVar_results)

## gnomAD, variant af and depth filtering 
clinVar_results <- clinVar_results %>% 
   mutate(variant_depth = if_else( as.integer( str_match(INFO, "DP\\=(\\d+)")[, 2])  > filter_variant_depth, "PASS","FAIL")) %>% 
   mutate(gnomad_af     = if_else( as.numeric( str_match(INFO, "gnomad_3_1_1_AF_non_cancer\\=(0\\.\\d+)\\;")[,2])  > filter_variant_af, "PASS","FAIL")) %>% 
   mutate(variant_af    = if_else(as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3]) / ( (as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,2]) ) + as.integer(str_match(Sample, ":(\\d+)\\,(\\d+)") [,3] )) > filter_variant_af, "PASS", "FAIL"))


## add column "vcf_id" to clinVar results in order to cross-reference with intervar and autopvs1 table
clinvar_results <- clinVar_results %>%
  mutate(c_id = str_match(INFO, "MedGen:(C\\d+);")[, 2]) %>% 
  mutate(vcf_id= str_remove_all(paste (CHROM,"-",START,"-",REF,"-",ALT), " "),
         ## add star annotations to clinVar results table based on filters // ## default version
         Stars = ifelse(grepl('CLNREVSTAT\\=criteria_provided,_single_submitter', INFO), "1",
                 ifelse(grepl('CLNREVSTAT\\=criteria_provided,_multiple_submitters', INFO), "2",
                 ifelse(grepl('CLNREVSTAT\\=reviewed_by_expert_panel', INFO), "3",
                 ifelse(grepl('CLNREVSTAT\\=practice_guideline', INFO), "4",
                 ifelse(grepl('CLNREVSTAT\\=criteria_provided,_conflicting_interpretations', INFO), "Needs Resolving", "0")
                                      )))),
          ## extract the calls and put in own column
         final_call = str_match(INFO, "CLNSIG\\=(\\w+)\\;")[, 2])

## retrieve and store clinVar input file into table data.table::fread()
input_submissions_file_path = file.path(input_dir, "submission_summary.txt")
 
submission_info_tab  <-  vroom(input_submissions_file_path, comment = "#",delim="\t", 
                               col_names = c("VariationID","ClinicalSignificance","DateLastEvaluated","Description","SubmittedPhenotypeInfo","ReportedPhenotypeInfo",
                                             "ReviewStatus","CollectionMethod","OriginCounts","Submitter","SCV","SubmittedGeneSymbol","ExplanationOfInterpretation"), 
                               show_col_types = FALSE)

submission_info_tab <- submission_info_tab %>% mutate(c_id = str_match(ReportedPhenotypeInfo, "(C\\d+):")[, 2])
       
## filter only those variants that need consensus call and find final call in submission table
entries_for_cc <- clinvar_results %>%
  filter(Stars == "Needs Resolving", na.rm = TRUE) 

entries_for_cc <- entries_for_cc %>% inner_join(submission_info_tab, by="c_id") %>% mutate(final_call=submission_info_tab$ClinicalSignificance)

## one Star cases that are “criteria_provided,_single_submitter” that do NOT have the B, LB, P, LP call must also go to intervar
additional_intervar_cases <- clinvar_results %>% filter(Stars == "1", final_call!="Benign",final_call!="Pathogenic", final_call != "Likely_benign",final_call!="Likeley_pathogenic")

## filter only those variant entries that need an InterVar run (No Star) and add the additional intervar cases from above
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

## add column "vcf_id" to intervar results in order to match with autopvs1 and clinvar table
intervar_results <- intervar_results %>%
                    mutate(vcf_id = str_remove_all(paste (Chr,"-",Start,"-",Ref,"-",Alt), " "))

## join all three tables together based on variant id that need intervar run
combined_tab_for_intervar <- autopvs1_results %>%
  inner_join(intervar_results, by="vcf_id") %>%
  filter(vcf_id %in% entries_for_intervar$vcf_id) %>%

  ## add column for individual scores that will be re-calculated if we need to adjust using autoPVS1 result
  mutate(evidencePVS1 = str_match(`InterVar: InterVar and Evidence`, "PVS1\\=(\\d+)\\s")[, 2]) %>%
  mutate(evidenceBA1 = str_match(`InterVar: InterVar and Evidence`, "BA1\\=(\\d+)\\s")[, 2]) %>%
  mutate( evidencePS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidencePM = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPM\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidencePP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sPP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidenceBS = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBS\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%
  mutate( evidenceBP = map_dbl(str_match(`InterVar: InterVar and Evidence`, "\\sBP\\=\\[([^]]+)\\]")[, 2], function(x) sum(as.integer(unlist(str_split(x, ",")))))     ) %>%

  ## if PVS1=0, take intervar call as final call
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
  mutate(evidencePP   = if_else((criterion == "IC4") & evidencePVS1 == 1, as.numeric(evidencePP)+1, as.double(evidencePP))) %>%

  #if criterion is na then PVS1 = 0;
  mutate(evidencePVS1 = if_else( (criterion == "na") & evidencePVS1 == 1, 0, as.double(evidencePVS1)))

## add new call based on new re-calculated scores (New ClinSig)
# Pathogenic - Criteria 1
# (i) 1 Very strong (PVS1) AND
# (a) ≥1 Strong (PS1–PS4) OR
# (b) ≥2 Moderate (PM1–PM6) OR
# (c) 1 Moderate (PM1–PM6) and 1 supporting (PP1–PP5) OR
# (d) ≥2 Supporting (PP1–PP5)
combined_tab_for_intervar <- combined_tab_for_intervar %>%
  mutate(final_call = if_else( (evidencePVS1 == 1) &
                               (evidencePS   >= 1) |
                               (evidencePM   >=2 ) |
                               (evidencePM   >= 1 & evidencePP ==1) |
                               (evidencePP   >=2 ), "Pathogenic", "")) %>%

# Pathogenic - Criteria 2
# (ii) ≥2 Strong (PS1–PS4) OR
  mutate(final_call = if_else( (evidencePS >= 2), "Pathogenic", "Uncertain significance"))  %>%

# Pathogenic - Criteria 3
# (iii) 1 Strong (PS1–PS4) AND
# (a)≥3 Moderate (PM1–PM6) OR
# (b)2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
# (c)1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5)
  mutate(final_call = if_else( (evidencePS == 1) &
                                 (evidencePM   >= 3) |
                                 (evidencePM   ==2 & evidencePP >=2 ) |
                                 (evidencePM == 1 & evidencePP >=4 ), "Pathogenic", "Uncertain significance")) %>%


# Likely pathogenic
# (i) 1 Very strong (PVS1) AND 1 moderate (PM1– PM6) OR
# (ii) 1 Strong (PS1–PS4) AND 1–2 moderate (PM1–PM6) OR
# (iii) 1 Strong (PS1–PS4) AND ≥2 supporting (PP1–PP5) OR
# (iv) ≥3 Moderate (PM1–PM6) OR
# (v) 2 Moderate (PM1–PM6) AND ≥2 supporting (PP1–PP5) OR
# (vi) 1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5)
  mutate(final_call = if_else( (evidencePVS1 == 1 & evidencePM == 1) |
                             (evidencePS==1 & evidencePM >= 1) |
                             (evidencePS==1 & evidencePP >=2 ) |
                             (evidencePM >= 3) |
                             (evidencePM ==2 & evidencePP>=2 ) |
                             (evidencePM == 1 & evidencePP>=4), "Likely Pathogenic", "Uncertain significance")) %>%

# Benign
# (i) 1 Stand-alone (BA1) OR
# (ii) ≥2 Strong (BS1–BS4)
  mutate(final_call = if_else( (evidenceBA1 == 1) |
                             (evidenceBS   >= 2), "Pathogenic", "Uncertain significance")) %>%

  # Likely Benign
# (i) 1 Strong (BS1–BS4) and 1 supporting (BP1– BP7) OR
# (ii) ≥2 Supporting (BP1–BP7)
  mutate(final_call = if_else( (evidenceBS == 1 & evidenceBP == 1) |
                             (evidenceBP   >= 2), "Pathogenic", "Uncertain significance"))

# Uncertain significance
# (i) non of the criteria were met.
# (ii) Benign and pathogenic are contradictory.

## merge tables together (clinvar and intervar) and write to file
master_tab <- full_join(clinvar_results,intervar_results, by= "vcf_id" )
write.table(master_tab, output_tab_file, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, quote = FALSE, col.names = TRUE)

## abridged version
results_tab_abridged <- master_tab %>% select(vcf_id,Gene.ensGene,gnomad_af, variant_depth, variant_af, final_call)
write.table(results_tab_abridged, output_tab_abr_file, append = FALSE, sep = "\t", dec = ".",
           row.names = FALSE, quote = FALSE, col.names = TRUE)