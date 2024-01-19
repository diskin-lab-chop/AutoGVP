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
    --workflow*|-w*)
      if [[ "$1" != *=* ]]; then shift; fi # Value is next arg if no `=`
      workflow="${1#*=}"
      ;;
    --vcf*|-v*)
      if [[ "$1" != *=* ]]; then shift; fi
      vcf_file="${1#*=}"
      ;;
    --filter_criteria*|-f*)
      if [[ "$1" != *=* ]]; then shift; fi
      filtering_criteria="${1#*=}"
      ;;
    --clinvar*|-c*)
      if [[ "$1" != *=* ]]; then shift; fi
      clinvar_file="${1#*=}"
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
    --selected_clinvar_submissions*)
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
    --help|-h)
        echo "Usage: $0 [-w/--workflow] [-v/--vcf <*.vcf*>] [-f/--filter_criteria=<criteria>] [-c/--clinvar <*.vcf>] [-i/--intervar <*.txt.intervar>] [-a/--autopvs1 <*autopvs1*>] [-m/--multianno <*multianno*>] [-o/--out <output>]"
        echo ""
        echo "Options:"
        echo "  -w/--workflow                     workflow"
        echo "  -v/--vcf                          VCF file"
        echo "  -f/--filter_criteria              VCF filtering criteria"
        echo "  -c/--clinvar                      clinvar file"
        echo "  -i/--intervar                     intervar results file"
        echo "  -a/--autopvs1                     autopvs1 results file"
        echo "  -m/--multianno                    multianno file"
        echo "  -O/--outdir                       output directory"
        echo "  -o/--out                          output prefix"
        echo "  --selected_clinvar_submissions    clinvar variant file with conflicts resolved"
        echo "  --variant_summary                 ClinVar variant summary file"
        echo "  --submission_summary              ClinVar submission summary file"
        echo "  --conceptIDs                      list of conceptIDs to prioritize submissions for clinvar variant conflict resolution. Will be ignored if selected_clinvar_submissions is provided"
        echo "  --conflict_res                    how to resolve conflicts associated with conceptIDs. Will be ignored if selected_clinvar_submissions is provided or if conceptIDs are not provided"
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

# Run appropriate Rscript based on workflow source (Cavatica vs. custom)
if [[ "$workflow" = "cavatica" ]] ; then

  if [[ -f $clinvar_file ]] ; then

    # Run AutoGVP from Cavatica workflow with provided clinvar vcf
    Rscript $BASEDIR/scripts/02-annotate_variants_CAVATICA_input.R --vcf $autogvp_input \
    --clinvar $clinvar_file \
    --multianno $multianno_input \
    --intervar $intervar_input \
    --autopvs1 $autopvs1_input \
    --output $out_file \
    --outdir $out_dir \
    --variant_summary $selected_submissions
    
    else
    
      # Run AutoGVP from Cavatica workflow with clinvar annotation in sample vcf
    Rscript $BASEDIR/scripts/02-annotate_variants_CAVATICA_input.R --vcf $autogvp_input \
    --multianno $multianno_input \
    --intervar $intervar_input \
    --autopvs1 $autopvs1_input \
    --output $out_file \
    --outdir $out_dir \
    --variant_summary $selected_submissions
    
  fi
  
  autogvp_output=${out_dir}/${out_file}".cavatica_input.annotations_report.abridged.tsv"

  else

  # Run AutoGVP from custom workflow
  Rscript $BASEDIR/scripts/02-annotate_variants_custom_input.R --vcf $autogvp_input \
  --clinvar $clinvar_file \
  --multianno $multianno_input \
  --intervar $intervar_input \
  --autopvs1 $autopvs1_input \
  --output $out_file \
  --outdir $out_dir \
  --variant_summary $selected_submissions \
  
  autogvp_output=${out_dir}/${out_file}".custom_input.annotations_report.abridged.tsv"

fi


# Parse vcf file so that info field values are in distinct columns
echo "Parsing VCF..."

bash $BASEDIR/scripts/03-parse_vcf.sh $autogvp_input

# Define parsed vcf and autogvp output file variables
vcf_parsed_file=${autogvp_input%.vcf*}."parsed.tsv"


# Filter VCF VEP gene/transcript annotations and merge data with AutoGVP output
echo "Filtering VEP annotations and creating final output..."

Rscript $BASEDIR/scripts/04-filter_gene_annotations.R --vcf $vcf_parsed_file --autogvp $autogvp_output --output $out_file --outdir $out_dir --colnames $BASEDIR/data/output_colnames.tsv

# Remove intermediate files
rm $autogvp_input $vcf_parsed_file $autogvp_output $out_dir/$out_file.filtered_csq_subfields.tsv $out_dir/${out_file}_multianno_filtered.txt $out_dir/${out_file}_autopvs1_filtered.tsv $out_dir/${out_file}_intervar_filtered.txt
