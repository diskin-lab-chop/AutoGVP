#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

echo $script_directory
input_file="$script_directory"/""

#-vcf input/clinvar_dummy.nih.norm.annot.chr17.vcf --intervar clinvar_input_wf_test.hg38_multianno.chr17.txt.intervar --autopvs1 clinvar_dummy.nih.autopvs1.chr17.tsv
# usage() { echo "Usage: $0 [-v <*.vcf>][-i <*.txt.intervar>] [-a <*autopvs1.tsv>]
#   [-w <cavatica || user] [-g <gnomad_var>] [-f genomAD_AF_filter <default: 0.001] [-v variant_depth_filter default: 15 ][-r variant_AF <default: 0.2>]" 1>&2; exit 1; }

## usage message to print for options
usage() {
  echo "Usage: $0 [-v <*.vcf>][-i <*.txt.intervar>] [-a <*autopvs1.tsv>]
  [-w <cavatica || user] [-g <gnomad_var>] [-f genomAD_AF_filter <default: 0.001] [-v variant_depth_filter default: 15 ][-r variant_AF <default: 0.2>]"
  echo ""
  echo "Options:"
  echo "  -v    vcf file"
  echo "  -i    intervar results file"
  echo "  -a    autopvs1 results file"
  echo "  -w    workflow type, must be either "cavatica" or "user""
  echo "  -g    gnomAD variable, default: gnomad_3_1_1_AF_non_cancer"
  echo "  -f    gnomAD allele frequency filter, default: 0.001)"
  echo "  -v    variant depth filter (default: 15)"
  echo "  -r    variant_AF <variant allele frequency, default: 0.2"
  echo "  -h    Display usage information."
  1>&2; exit 1; }

## default values for options
clinvar_version="clinvar_20211225"
genomAD_AF_filter=0.001
variant_depth_filter=15
variant_AF=.2
workflow_type="user"

## if cavatica worklflow save gnomad variable as "gnomad_3_1_1_AF_non_cancer"
if [ "$workflow_type" == 'cavatica' ]
then
  gnomad_var="gnomad_3_1_1_AF_non_cancer"
fi

while getopts ":v:i:a:w:g:f:v:r:h" arg; do
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
        w) #workflow type
          workflow_type="$OPTARG"
          ;;
        g) #gnomad variable/column
          gnomad_var="$OPTARG"
          ;;
        f) #genomAD_AF_filter
          genomAD_AF_filter="$OPTARG"
          ;;
        v) #variant_depth_filter type
          variant_depth_filter="$OPTARG"
          ;;
        r) #variant_AF type
          variant_AF="$OPTARG"
          ;;
        h | *) # Display help.
          usage
          exit 0
          ;;
    esac
done


echo "vcf = $vcf_file"
echo "intervar = $intervar_file"
echo "autopvs1 = $autopvs1_file"
echo "workflow_type = $workflow_type"
echo "gnomad_var= $gnomad_var"
echo "genomAD_AF_filter= $genomAD_AF_filter"
echo "variant_depth_filter = $variant_depth_filter"
echo "variant_AF = $variant_AF"
echo "workflow_type = $workflow_type"

if [[ "$workflow_type" == 'user' && -z ${gnomad_var} ]]
then
  echo "ERROR: if workflow type is non-cavatica, must provide gnomAD_var (ie. 'gnomad_3_1_1_AF_non_cancer') ";
  exit 1;
fi

## if workflow is cavatica run Rscript with this cmd
if [ "$workflow_type" == 'cavatica' ]
then
  echo "Rscript 01-annotate_variants.R --vcf $vcf_file --intervar intervar_file --autopvs1 autopvs1_file --clinvar $clinvar_version --gnomad_variable $gnomad_var --gnomad_af $genomAD_AF_filter --variant_depth $variant_depth_filter --variant_af variant_AF"
  pwd
fi

##default file paths for histology and rmats output
#hist_file="../../data/v19_plus_20210311_pnoc_rna.tsv"
#rmats_file="../../data/merge_rMATS_splicing.SE.single.tsv"

## run script
#Rscript 01-annotate_variants.R --vcf input/clinvar_dummy.nih.norm.annot.chr17.vcf --intervar clinvar_input_wf_test.hg38_multianno.chr17.txt.intervar --autopvs1 clinvar_dummy.nih.autopvs1.chr17.tsv
