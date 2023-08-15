#!/usr/bin/env bash

# Define variables
vcf_file=$1
out_file=$2
out_dir=$3
exp_args=("$@")


# Define filtered vcf output file 
vcf_filtered_file=${out_file}."filtered.vcf"

# Print input vcf file name
echo "vcf file: $vcf_file ";

# define default filters
cmd="bcftools view -f 'PASS,.' $vcf_file"

## loop through args for other user-defined filters
for i in "${exp_args[@]:3}"; do
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
