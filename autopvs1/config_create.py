"""
A silly helper script to createa config file at run time based on where input references directry is.
Ref directory should be in PWD, config file will be created in PWD
This tool is most useful for use in automated/non-interactive envrinoments
"""
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', help='Name of dir containing reference files. Should be a subdir in the cwd')
args = parser.parse_args()

data_dir = args.data_dir
cwd = os.getcwd()
"""
Can be made more flexible in the future.
Config file format will basically be this, with cwd and data dir paths bein g used with fixed file names:

[DEFAULT]
pvs1levels = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/PVS1.level
gene_alias = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/hgnc.symbol.previous.tsv
gene_trans = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/clinvar_trans_stats.tsv

[HG38]
genome = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/hg38.fa
transcript = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/ncbiRefSeq_hg38.gpe
domain = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/functional_domains_hg38.bed
hotspot = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/mutational_hotspots_hg38.bed
curated_region = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/expert_curated_domains_hg38.bed
exon_lof_popmax = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/exon_lof_popmax_hg38.bed
pathogenic_site = /home/ubuntu/tools/pathogenicity-assessment/autopvs1/data/clinvar_pathogenic_GRCh38.vcf
"""
# find fasta file used
contents = os.listdir(data_dir)
fasta = ""
for fname in contents:
    parts = fname.split('.')
    if parts[-1] == 'fasta' or parts[-1] == 'fa':
        fasta = fname
        break
if len(fasta) == 0:
    sys.stderr.write("ERROR: Fasta file not found. found these files: " + "\n".join(contents))
    exit(1)
# output config file....hardcoded. woof.
config_file = open('config.ini', 'w')
config_file.write(
    "[DEFAULT]\n"
    "pvs1levels = " + cwd + "/" + data_dir + "/PVS1.level\n"
    "gene_alias = " + cwd + "/" + data_dir + "/hgnc.symbol.previous.tsv\n"
    "gene_trans = " + cwd + "/" + data_dir + "/clinvar_trans_stats.tsv\n\n"
    "[HG38]\n"
    "genome = " + cwd + "/" + data_dir + "/" + fasta + "\n"
    "transcript = " + cwd + "/" + data_dir + "/ncbiRefSeq_hg38.gpe\n"
    "domain = " + cwd + "/" + data_dir + "/functional_domains_hg38.bed\n"
    "hotspot = " + cwd + "/" + data_dir + "/mutational_hotspots_hg38.bed\n"
    "curated_region = " + cwd + "/" + data_dir + "/expert_curated_domains_hg38.bed\n"
    "exon_lof_popmax = " + cwd + "/" + data_dir + "/exon_lof_popmax_hg38.bed\n"
    "pathogenic_site = " + cwd + "/" + data_dir + "/clinvar_pathogenic_GRCh38.vcf\n"
     )

config_file.close()