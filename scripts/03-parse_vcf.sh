#!/usr/bin/env bash
vcf_file=$1

# define parsed output file and csq field file variables
vcf_parsed_file=${vcf_file%.vcf*}."parsed.tsv"
csq_fields=${vcf_file%.vcf*}_"csq_subfields.tsv"

# Extract list of subfields in INFO column; if multiple, these will be separated by a semicolon
egrep -v "^#" $vcf_file | awk '{ n=split($8, tmp, /=[^;]*;*|;/); for(i=1; i<n; i++) print tmp[i] }' | sort -u > ${vcf_file%.vcf*}_subfields_all.tsv

# If only one INFO subfield, it will precede "="
egrep -v "^#" $vcf_file | awk '{print $8}' | sed 's/=.*//' | sort -u > ${vcf_file%.vcf*}_subfields_first.tsv

# Merge extracted subfields and assign to variable
cat ${vcf_file%.vcf*}_subfields_all.tsv ${vcf_file%.vcf*}_subfields_first.tsv | sort -u > ${vcf_file%.vcf*}_subfields.tsv
subfields=$(cat ${vcf_file%.vcf*}_subfields.tsv)

# Call function to join subfields
function join_by {
  local d=${1-} f=${2-}
  if shift 2; then
    printf %s "$f" "${@/#/$d}"
  fi
}

# Join subfields by '\t%' for bcftools query
subfield_list=$(join_by '\t%' $subfields)

# Add standard vcf columns to list to parse
all_columns="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%${subfield_list}\n"

# run vcftools query to parse all columns and info subfields, and include header
bcftools query -H -f $all_columns $vcf_file > $vcf_parsed_file

# Extract VEP CSQ field names and save to tsv
bcftools +split-vep $vcf_file -l | awk -F '\t' '{print $2}' > $csq_fields

# remove intermediate files
rm ${vcf_file%.vcf*}_subfields.tsv ${vcf_file%.vcf*}_subfields_first.tsv ${vcf_file%.vcf*}_subfields_all.tsv
