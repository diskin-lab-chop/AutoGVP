#!/usr/bin/env bash

set -e

# export bcftools plugin environmental variable
export BCFTOOLS_PLUGINS=/rocker-build/bcftools-1.17/plugins

# Define root directory of repo
BASEDIR="$(dirname "${BASH_SOURCE[0]}")"
echo "$BASEDIR"

## default files
variant_summary_file="$BASEDIR/data/ClinVar-selected-submissions.tsv"

# define parameter variables 
while [ $# -gt 0 ]; do
  case "$1" in
    --vcf*|-v*)
      if [[ "$1" != *=* ]]; then shift; fi
      vcf_file="${1#*=}"
      ;;
    --filter_criteria*|-f*)
      if [[ "$1" != *=* ]]; then shift; fi
      filtering_criteria="${1#*=}"
      ;;
    --intervar*|-i*)
      if [[ "$1" != *=* ]]; then shift; fi
      intervar_file="${1#*=}"
      ;;
    --multianno*|-m*)
      if [[ "$1" != *=* ]]; then shift; fi
      multianno_file="${1#*=}"
      ;;
    --autopvs1*|-a*)
      if [[ "$1" != *=* ]]; then shift; fi
      autopvs1_file="${1#*=}"
      ;;
    --outdir*|-O*)
      if [[ "$1" != *=* ]]; then shift; fi
      out_dir="${1#*=}"
      ;;
    --out*|-o*)
      if [[ "$1" != *=* ]]; then shift; fi
      out_file="${1#*=}"
      ;;
    --sample_id*|-s*)
      if [[ "$1" != *=* ]]; then shift; fi
      sample_id="${1#*=}"
      ;;
    --selected_submissions*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      selected_submissions="${1#*=}"
      ;;
    --variant_summary*)
      if [[ "$1" != *=* ]]; then shift; fi
      variant_summary="${1#*=}"
      ;;
    --submission_summary*)
      if [[ "$1" != *=* ]]; then shift; fi
      submission_summary="${1#*=}"
      ;;
    --conceptIDs*)
      if [[ "$1" != *=* ]]; then shift; fi
      conceptIDs="${1#*=}"
      ;;
    --conflict_res*)
      if [[ "$1" != *=* ]]; then shift; fi
      conflict_res="${1#*=}"
      ;;
    --custom_output_cols*)
      if [[ "$1" != *=* ]]; then shift; fi
      custom_output_cols="${1#*=}"
      ;;
    --help|-h)
        echo "Usage: $0 [-v/--vcf <*.vcf*>] [-f/--filter_criteria=<criteria>] [-i/--intervar <*.txt.intervar>] [-a/--autopvs1 <*autopvs1*>] [-m/--multianno <*multianno*>] [-c/--selected_clinvar_submissions <refs/ClinVar-selected-submissions.tsv>] [-o/--out <output>]"
        echo ""
        echo "Options:"
        echo "  -v/--vcf                          VCF file"
        echo "  -f/--filter_criteria              VCF filtering criteria"
        echo "  -i/--intervar                     Intervar results file"
        echo "  -a/--autopvs1                     Autopvs1 results file"
        echo "  -m/--multianno                    ANNOVAR file"
        echo "  -O/--outdir                       output directory"
        echo "  -o/--out                          output prefix"
        echo "  -s/--sample_id                    sample ID to be added to the output file"
        echo "  -c/--selected_submissions         ClinVar variant file with conflicts resolved"
        echo "  --variant_summary                 ClinVar variant summary file"
        echo "  --submission_summary              ClinVar submission summary file"
        echo "  --conceptIDs                      list of conceptIDs to prioritize submissions for clinvar variant conflict resolution. Will be ignored if selected_clinvar_submissions is provided"
        echo "  --conflict_res                    how to resolve conflicts associated with conceptIDs. Will be ignored if selected_clinvar_submissions is provided or if conceptIDs are not provided"
        echo "  --custom_output_cols              optional; text file of user-defined column names from VCF info fields or other input file to be included in AutoGVP output files. Must contain three columns named 'Column_name', 'Rename' (i.e., what to rename colum in final output), and 'Abridged' (T or F indicating if column should be included in abridged output)"
        echo "  -h/--help                         Display usage information"
            exit 0
      ;;
    *)
      >&2 printf "Error: Invalid argument\n"
      exit 1
      ;;
  esac
  shift
done


# If selected ClinVar submissions files not provided, then run select-ClinVar-submissions.R
if [[ ! -e $selected_submissions ]]; then
    
    echo "select ClinVar submission file not specified. Running select-ClinVar-submissions Rscript..."
    
  # if variant_summary or submission_summary args not provided, check if files exist in data/
#  if [[ $variant_summary -eq 0  || $submission_summary -eq 0 ]]
  if [[ ! -e $variant_summary || ! -e $submission_summary ]]
  then
      
      echo "variant summary and/or submission_summary file(s) not specified. Checking if files exist in data/..."
      
      variant_summary=$(find data/ -type f -name "variant_summary*")
      submission_summary=$(find data/ -type f -name "submission_summary*")
      
      # if no files found matching pattern, download latest versions from ClinVar
      # if [[ -z "$variant_summary" || -z "$submission_summary" ]]
      if [[ ! -e "$variant_summary" || ! -e "$submission_summary" ]]
      then
        
        echo "variant_summary and/or submission_summary files not found. Downloading latest versions from ClinVar..."
        
        bash $BASEDIR/scripts/download_db_files.sh
        
        variant_summary=$(find data/ -type f -name "variant_summary*")
        submission_summary=$(find data/ -type f -name "submission_summary*") ; 
        
      else 
        
        echo "variant and submission summary files found. Running select-ClinVar-submissions Rscript..."
        
      fi
      
    fi
      
        if [[ ! -e $conceptIDs ]] ; then
        
        echo "resolving ClinVar conflicts using default parameters..."
      
        Rscript $BASEDIR/scripts/select-clinVar-submissions.R --variant_summary $variant_summary --submission_summary $submission_summary --outdir $out_dir
        
      fi
      
      if [[ -f $conceptIDs && -z $conflict_res ]] ; then 
      
        echo "resolving ClinVar conflicts with provided concept IDs and taking latest date evaluated call..."
      
        Rscript $BASEDIR/scripts/select-clinVar-submissions.R --variant_summary $variant_summary --submission_summary $submission_summary --outdir $out_dir --conceptID_list $conceptIDs --conflict_res "latest"
        
      fi
      
      if [[ -f $conceptIDs && -n $conflict_res ]] ; then
      
        echo "resolving ClinVar conflicts with provided concept IDs and specified conflict resolution..."
      
        Rscript $BASEDIR/scripts/select-clinVar-submissions.R --variant_summary $variant_summary --submission_summary $submission_summary --outdir $out_dir --conceptID_list $conceptIDs --conflict_res $conflict_res
      
      fi
      
      selected_submissions="$BASEDIR/results/ClinVar-selected-submissions.tsv"
        
fi


# Filter VCF file; by default the function performs filtering based on FILTER column, with other criteria specified by user
echo "Filtering VCF..."

vcf_filtered_file=${out_file}."filtered.vcf"
bash $BASEDIR/scripts/01-filter_vcf.sh $vcf_file $multianno_file $autopvs1_file $intervar_file $out_file $out_dir $filtering_criteria

autogvp_input=$out_dir/$vcf_filtered_file
multianno_input=$out_dir/${out_file}_multianno_filtered.txt
autopvs1_input=$out_dir/${out_file}_autopvs1_filtered.tsv
intervar_input=$out_dir/${out_file}_intervar_filtered.txt

# Run AutoGVP 
echo "Running AutoGVP..."

# Run AutoGVP variant annotation
Rscript $BASEDIR/scripts/02-annotate_variants.R --vcf $autogvp_input \
  --clinvar $selected_submissions \
  --multianno $multianno_input \
  --intervar $intervar_input \
  --autopvs1 $autopvs1_input \
  --output $out_file \
  --outdir $out_dir \
  --sample_id $sample_id
  
autogvp_output=${out_dir}/${out_file}".annotations_report.abridged.tsv"


# Parse vcf file so that info field values are in distinct columns
echo "Parsing VCF..."

bash $BASEDIR/scripts/03-parse_vcf.sh $autogvp_input

# Define parsed vcf and autogvp output file variables
vcf_parsed_file=${autogvp_input%.vcf*}."parsed.tsv"


# Filter VCF VEP gene/transcript annotations and merge data with AutoGVP output
echo "Filtering VEP annotations and creating final output..."

  # check if custom colnames file is provided
  if [[ -f $custom_output_cols ]] ; then

  Rscript $BASEDIR/scripts/04-filter_gene_annotations.R --vcf $vcf_parsed_file \
  --autogvp $autogvp_output \
  --output $out_file \
  --outdir $out_dir \
  --default_colnames $BASEDIR/data/output_colnames_default.tsv \
  --custom_colnames $custom_output_cols

  else 

  Rscript $BASEDIR/scripts/04-filter_gene_annotations.R --vcf $vcf_parsed_file \
  --autogvp $autogvp_output \
  --output $out_file \
  --outdir $out_dir \
  --default_colnames $BASEDIR/data/output_colnames_default.tsv

  fi

# Remove intermediate files
rm $autogvp_input $vcf_parsed_file $autogvp_output $out_dir/$out_file.filtered_csq_subfields.tsv $out_dir/${out_file}_multianno_filtered.txt $out_dir/${out_file}_autopvs1_filtered.tsv $out_dir/${out_file}_intervar_filtered.txt
