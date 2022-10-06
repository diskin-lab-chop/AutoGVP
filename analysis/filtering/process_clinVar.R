library("dplyr")
library("tidyverse")

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

## set up directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses")

input_dir   <- file.path(root_dir, "input")

## retrieve and store clinVar input file into table
input_clinVar_file <- (paste0(input_dir, "/", "clinvar_dummy.nih.norm.annot.vcf.gz"))
clinVar_results  <-  read.delim(input_clinVar_file, skip = 127, sep = "\t", header=TRUE) 

#test_file <- (paste0(input_dir, "/", "clinvar_dummy.nih.test.tsv"))
#clinVar_test <-  read.delim(test_file, skip = 127, sep = "\t", header=TRUE) 

## add star annotations to clinVar results table based on filters
clinVar_results_stars <- mutate(clinVar_results, 
                                Stars = ifelse(grepl('CLNREVSTAT\\=criteria_provided,_single_submitter', INFO), "One",
                                         ifelse(grepl('CLNREVSTAT\\=criteria_provided,_multiple_submitters', INFO), "Two",
                                          ifelse(grepl('CLNREVSTAT\\=reviewed_by_expert_panel', INFO), "Three",
                                           ifelse(grepl('CLNREVSTAT\\=practice_guideline', INFO), "Four",
                                            ifelse(grepl('CLNREVSTAT\\=criteria_provided,_conflicting_interpretations', INFO), "Needs Resolving", "No Star")
                                            ))))) %>% 
  
                          ## extract the calls and put in own column
                          mutate(Call = str_match(INFO, "CLNSIG\\=(\\w+)\\;")[, 2])

 

## filter only those variants that need an InterVar run (No Star)
run_interVar <- filter(clinVar_results_stars, Stars == "No Star", na.rm = TRUE)

## filter only those variants that need consensus call
run_cc <- filter(clinVar_results_stars, Stars == "Needs Resolving", na.rm = TRUE)


