#!/usr/bin/env bash
vcf_file=$1

## name output file
vcf_parsed_file=${vcf_file%.vcf*}."parsed.tsv"

## Extract list of subfields in INFO column
egrep -v "^#" $vcf_file | awk '{ n=split($8, tmp, /=[^;]*;/); for(i=1; i<n; i++) print tmp[i] }' | sort -u > subfields.tsv

subfields=$(cat subfields.tsv)

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
rm subfields.tsv
