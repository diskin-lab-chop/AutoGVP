#!/usr/bin/env bash

# Define variables
vcf_file=$1
multianno_file=$2
autopvs1_file=$3
intervar_file=$4
out_file=$5
out_dir=$6
exp_args=("$@")

# Define filtered vcf output file 
vcf_filtered_file=${out_file}."filtered.vcf"

# Print input vcf file name
echo "vcf file: $vcf_file ";

# define default filters
cmd="bcftools view -f 'PASS,.' $vcf_file"

## loop through args for other user-defined filters
for i in "${exp_args[@]:6}"; do
  #echo filtering for... $i
  cmd+=" | bcftools filter -i '$i'"
#  cmd+=" $vcf_file "
  #echo $cmd
done

# define full bcftools filter command
cmd+=" > $out_dir/$vcf_filtered_file"

# print command
echo "cmd: " $cmd

# execute command
eval "$cmd"

# extract variant keys (chrom-pos-ref-alt, without "chr") from the filtered vcf; these match
# the `vcf_id` built in 02-annotate_variants.R. Matching on the full key (rather than on
# position alone) avoids retaining unrelated rows that merely contain the same number.
vcf_keys_file=$out_dir/${out_file}_vcf_filtered_keys.tsv
bcftools query -f "%CHROM-%POS-%REF-%ALT\n" $out_dir/$vcf_filtered_file | sed 's/^chr//' > $vcf_keys_file

# filter multianno file for filtered vcf variants; the vcf key is Chr-Otherinfo5-Otherinfo7-Otherinfo8
echo "Filtering multianno file..."

multianno_filtered_file=${out_file}_multianno_filtered.txt
gzip -cdf $multianno_file | awk -F'\t' '
  FILENAME == ARGV[1] { keys[$1]; next }
  FNR == 1 {
    for (i = 1; i <= NF; i++) { if ($i == "Otherinfo5") p = i; if ($i == "Otherinfo7") r = i; if ($i == "Otherinfo8") a = i }
    if (!p || !r || !a) { print "ERROR: Otherinfo5/7/8 columns not found in multianno file" > "/dev/stderr"; exit 1 }
    print; next
  }
  { c = $1; sub(/^chr/, "", c); if ((c "-" $p "-" $r "-" $a) in keys) print }
' $vcf_keys_file - > $out_dir/$multianno_filtered_file

# filter autopvs1 file for filtered vcf variants; vcf_id is the first column
echo "Filtering autopvs1 file..."

autopvs1_filtered_file=${out_file}_autopvs1_filtered.tsv
gzip -cdf $autopvs1_file | awk -F'\t' '
  FILENAME == ARGV[1] { keys[$1]; next }
  FNR == 1 { print; next }
  { k = $1; gsub(/ /, "", k); sub(/^chr/, "", k); if (k in keys) print }
' $vcf_keys_file - > $out_dir/$autopvs1_filtered_file

# filter intervar file for the variants retained in the multianno file (Chr, Start, Ref, Alt)
echo "Filtering intervar file..."

intervar_filtered_file=${out_file}_intervar_filtered.txt
awk -F'\t' '
  FILENAME == ARGV[1] { c = $1; sub(/^chr/, "", c); keys[c "-" $2 "-" $4 "-" $5]; next }
  FNR == 1 { print; next }
  { c = $1; sub(/^chr/, "", c); if ((c "-" $2 "-" $4 "-" $5) in keys) print }
' $out_dir/$multianno_filtered_file <(gzip -cdf $intervar_file) > $out_dir/$intervar_filtered_file

# remove key file
rm $vcf_keys_file
