#!/usr/bin/env bash


## default files
submission_summary_file="input/submission_summary.txt.gz"
variant_summary_file="input/variant_summary.txt.gz"


## usage message to print for options
usage() {
  echo "Usage: $0 [-v <*.vcf*>][-i <*.txt.intervar>] [-a <*autopvs1*>] [-n <*multianno*>] [-o <output>]"
  echo ""
  echo "Options:"
  echo "  -v    vcf file"
  echo "  -i    intervar results file"
  echo "  -a    autopvs1 results file"
  echo "  -m    multianno file"
  echo "  -c    clinvar db file"
  echo "  -o    output prefix"
  echo "  -h    Display usage information."
  1>&2; exit 1; }

while getopts ":v:i:a:m:o:c:h" arg; do
    case "$arg" in
        v) # vcf file
          vcf_file="$OPTARG"
          ;;
        i) # intervar file
          intervar_file="$OPTARG"
          ;;
        a) # autopvs1 file
          autopvs1_file="$OPTARG"
          ;;
        m) # multianno file
          multianno_file="$OPTARG"
          ;;
        o) # output prefix
          out_file="$OPTARG"
          ;;
        c) ## clinvar file
          clinvar_file="$OPTARG"
          ;;
        h | *) # Display help.
        usage
          exit 0
          ;;
    esac
done

vcf_filtered_file=${vcf_file%.vcf*}."filtered.vcf"


echo "bash filter_vcf.sh $vcf_file"
echo "Rscript 01-annotate_variants_custom_input.R --vcf $vcf_filtered_file --clinvar $clinvar_file --multianno $multianno_file --intervar $intervar_file --autopvs1 $autopvs1_file --output $out_file --variant_summary $variant_summary_file --submission_summary $submission_summary_file"
echo "bash parse_vcf.sh $vcf_filtered_file"

bash filter_vcf.sh $vcf_file
Rscript 01-annotate_variants_custom_input.R --vcf $vcf_filtered_file --clinvar $clinvar_file --multianno $multianno_file --intervar $intervar_file --autopvs1 $autopvs1_file --output $out_file --variant_summary $variant_summary_file --submission_summary $submission_summary_file
bash parse_vcf.sh $vcf_filtered_file
