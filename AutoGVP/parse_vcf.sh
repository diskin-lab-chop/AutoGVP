#!/usr/bin/env bash
vcf_file=$1

## name output file
vcf_parsed_file=${vcf_file%.vcf*}."parsed.vcf"

egrep -v "^#" $vcf_file | head -n 2 | awk '{ n=split($8, tmp, /=[^;]*;/); for(i=1; i<n; i++) print tmp[i] }' > subfields.tsv

subfields=$(cat subfields.tsv)

function join_by {
  local d=${1-} f=${2-}
  if shift 2; then
    printf %s "$f" "${@/#/$d}"
  fi
}

subfield_list=$(join_by '\t%' $subfields)

all_columns="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%${subfield_list}\n"

bcftools query -H -f $all_columns $vcf_file > $vcf_parsed_file