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

# extract vcf positions
bcftools query -f "%POS\n" $out_dir/$vcf_filtered_file > $out_dir/${out_file}_vcf_filtered_positions.tsv
echo "Chr" >> $out_dir/${out_file}_vcf_filtered_positions.tsv
echo "vcf_id" >> $out_dir/${out_file}_vcf_filtered_positions.tsv

# filter multianno file for filtered vcf positions
echo "Filtering multianno file..."

multianno_filtered_file=${out_file}_multianno_filtered.txt
zgrep -Fw -f $out_dir/${out_file}_vcf_filtered_positions.tsv $multianno_file > $out_dir/$multianno_filtered_file

# filter autopvs1 file for filtered vcf positions
echo "Filtering autopvs1 file..."

autopvs1_filtered_file=${out_file}_autopvs1_filtered.tsv
zgrep -Fw -f $out_dir/${out_file}_vcf_filtered_positions.tsv $autopvs1_file > $out_dir/$autopvs1_filtered_file

# extract multianno file positions to match with intervar file positions
cat $out_dir/$multianno_filtered_file | awk '{print $2}' > $out_dir/${out_file}_multianno_filtered_positions.tsv
echo "#Chr" >> $out_dir/${out_file}_multianno_filtered_positions.tsv

# filter intervar file for filtered multianno positions
echo "Filtering intervar file..."

intervar_filtered_file=${out_file}_intervar_filtered.txt
zgrep -Fw -f $out_dir/${out_file}_multianno_filtered_positions.tsv $intervar_file > $out_dir/$intervar_filtered_file

# remove position files
rm $out_dir/${out_file}_vcf_filtered_positions.tsv $out_dir/${out_file}_multianno_filtered_positions.tsv 
