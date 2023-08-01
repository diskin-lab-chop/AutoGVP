#!/usr/bin/env bash


## default files
variant_summary_file="input/ClinVar-selected-submissions.tsv"
# submission_summary_file="input/submission_summary.txt.gz"
# variant_summary_file="input/variant_summary.txt.gz"

if [[ ! -f $variant_summary_file ]] ; then
    echo "ERROR: ClinVar-selected-submissions.tsv file not found. Please run select-clinvar-submissions.R script before running AutoGVP"
fi


## usage message to print for options
# usage() {
#   echo "Usage: $0 [-v <*.vcf*>][-i <*.txt.intervar>] [-a <*autopvs1*>] [-n <*multianno*>] [-o <output>] [-f <filtering_criteria>]"
#   echo ""
#   echo "Options:"
#   echo "  -v    vcf file"
#   echo "  -i    intervar results file"
#   echo "  -a    autopvs1 results file"
#   echo "  -m    multianno file"
#   echo "  -o    output prefix"
#   echo "  -f    additional filtering criteria"
#   echo "  -h    Display usage information."
#   1>&2; exit 1; }

while getopts ":v:i:a:m:o:h" arg; do
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
        h | *) # Display help.
        # usage
        #   exit 0
          ;;
    esac
done

vcf_filtered_file=${vcf_file%.vcf*}."filtered.vcf"

# Filter VCF file for read depth, alt allele freq, pass filter and other specified criteria
bash filter_vcf.sh $vcf_file

# Run AutoGVP from Cavatica workflow
Rscript 01-annotate_variants_CAVATICA_input.R --vcf $vcf_filtered_file --multianno $multianno_file --intervar $intervar_file --autopvs1 $autopvs1_file --output $out_file --variant_summary $variant_summary_file

# Parse vcf file so that info field values are in distinct columns
bash parse_vcf.sh $vcf_filtered_file

# Define parsed vcf and autogvp output file variables
vcf_parsed_file=${vcf_filtered_file%.vcf*}."parsed.tsv"
autogvp_output=${out_file}".cavatica_input.annotations_report.abridged.tsv"

# Filter VCF gene/transcript annotations and merge data with AutoGVP output
Rscript 04-filter_gene_annotations.R --vcf $vcf_parsed_file --autogvp $autogvp_output --output $out_file

rm $vcf_filtered_file $vcf_parsed_file
