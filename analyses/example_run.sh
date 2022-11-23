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

## run script
Rscript 01-annotate_variants.R --vcf input/clinvar_dummy.nih.norm.annot.chr17.vcf --intervar clinvar_input_wf_test.hg38_multianno.chr17.txt.intervar --autopvs1 clinvar_dummy.nih.autopvs1.chr17.tsv
