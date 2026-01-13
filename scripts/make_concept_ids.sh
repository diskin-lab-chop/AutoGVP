#!/usr/bin/env bash

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
BASEDIR="$(pwd)"

# Get date string
DATE="$(date +'%Y%m%d')"

# Set file urls
disease_file="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/disease_names"
gene_file="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"
medgen_file="https://ftp.ncbi.nlm.nih.gov/pub/medgen/MedGenIDMappings.txt.gz"
cpg_file="https://raw.githubusercontent.com/diskin-lab-chop/pbta-germline-somatic/refs/heads/main/analyses/oncokb-annotation/input/cpg.txt"
mondo_file="https://purl.obolibrary.org/obo/mondo.obo"

# Make full disease concept ID list
curl -fsSL $disease_file | awk -F'\t' '($7 == "Disease" && $1 != "not specified" && $1 != "not provided") {print $3}' | sort -u \
  > $BASEDIR/../data/clinvar_all_disease_concept_ids_$DATE.txt

# Make cpg concept ID list
awk -F'\t' '
  NR==FNR { list[$1]=1; next }
  ($2 in list) {print $4}
' <(curl -fsSL "$cpg_file") <(curl -fsSL "$gene_file") | sort -u \
  > $BASEDIR/../data/clinvar_cpg_concept_ids_$DATE.txt

# Get cancer associated MONDO terms (under MONDO:0045024)
curl -fsSL https://purl.obolibrary.org/obo/mondo.obo |
awk '
  /^\[Term\]/ { id=""; n=0; next }
  /^id: MONDO:/ { id=$2; next }
  /^is_a: MONDO:/ { parent[++n]=$2; next }
  /^$/ {
    if (id) for (i=1; i<=n; i++) print id "\t" parent[i]
  }
' |
awk -v root="MONDO:0045024" '
  BEGIN { keep[root]=1 }
  {
    child=$1; parent=$2
    # record edges so we can iterate until no new nodes are added
    kids[parent] = kids[parent] SUBSEP child
  }
  END {
    changed=1
    while (changed) {
      changed=0
      for (p in keep) {
        n = split(kids[p], arr, SUBSEP)
        for (i=1; i<=n; i++) {
          c = arr[i]
          if (c != "" && !(c in keep)) { keep[c]=1; changed=1 }
        }
      }
    }
    for (id in keep) print id
  }
' | sort -u > $BASEDIR/../data/mondo_cancer_$DATE.txt

# Make cancer-associated concept ID list
awk -F'\t' '
  NR==FNR { list[$1]=1; next }
  ($4 in list && $7 == "Disease") {print $3}
' $BASEDIR/../data/mondo_cancer_$DATE.txt <(curl -fsSL "$disease_file") | sort -u \
  > $BASEDIR/../data/clinvar_cancer_concept_ids_$DATE.txt
