# Pathogenicity assessment development

## Prequisite
It is recommended to have first run the [Kids First Germline Annotation Workflow](https://github.com/kids-first/kf-germline-workflow/blob/master/docs/GERMLINE_SNV_ANNOT_README.md) first.

## Pathogenicity Preprocessing Workflow
This workflow uses the prerequisite input to run the InverVar workflow and autoPVS1 tool.
Recommended inputs:
 - `annovar_db`: Annovar Database with at minimum required resources to InterVar. Need to use [annovar download commands](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) to get the following:
     ```
            annovar_humandb_hg38_intervar/
        ├── hg38_AFR.sites.2015_08.txt
        ├── hg38_AFR.sites.2015_08.txt.idx
        ├── hg38_ALL.sites.2015_08.txt
        ├── hg38_ALL.sites.2015_08.txt.idx
        ├── hg38_AMR.sites.2015_08.txt
        ├── hg38_AMR.sites.2015_08.txt.idx
        ├── hg38_EAS.sites.2015_08.txt
        ├── hg38_EAS.sites.2015_08.txt.idx
        ├── hg38_EUR.sites.2015_08.txt
        ├── hg38_EUR.sites.2015_08.txt.idx
        ├── hg38_SAS.sites.2015_08.txt
        ├── hg38_SAS.sites.2015_08.txt.idx
        ├── hg38_avsnp147.txt
        ├── hg38_avsnp147.txt.idx
        ├── hg38_clinvar_20210501.txt
        ├── hg38_clinvar_20210501.txt.idx
        ├── hg38_dbnsfp42a.txt
        ├── hg38_dbnsfp42a.txt.idx
        ├── hg38_dbscsnv11.txt
        ├── hg38_dbscsnv11.txt.idx
        ├── hg38_ensGene.txt
        ├── hg38_ensGeneMrna.fa
        ├── hg38_esp6500siv2_all.txt
        ├── hg38_esp6500siv2_all.txt.idx
        ├── hg38_gnomad_genome.txt
        ├── hg38_gnomad_genome.txt.idx
        ├── hg38_kgXref.txt
        ├── hg38_knownGene.txt
        ├── hg38_knownGeneMrna.fa
        ├── hg38_refGene.txt
        ├── hg38_refGeneMrna.fa
        ├── hg38_refGeneVersion.txt
        ├── hg38_rmsk.txt
        └── hg38_seq
            ├── annovar_downdb.log
            └── hg38.fa
    ```
 - `intervar_db`: InterVar Database from git repo + mim_genes.txt
 - `autopvs1_db`: git repo files plus a user-provided fasta reference. For hg38, recommend:
    ```
    data
        ├── Homo_sapiens_assembly38.fasta
        ├── Homo_sapiens_assembly38.fasta.fai
        ├── PVS1.level
        ├── clinvar_pathogenic_GRCh38.vcf
        ├── clinvar_trans_stats.tsv
        ├── exon_lof_popmax_hg38.bed
        ├── expert_curated_domains_hg38.bed
        ├── functional_domains_hg38.bed
        ├── hgnc.symbol.previous.tsv
        ├── mutational_hotspots_hg38.bed
        └── ncbiRefSeq_hg38.gpe
    ```
 - `annovar_db_str`: Name of dir created when `annovar_db` tar ball in decompressed. Default: `annovar_humandb_hg38_intervar`
 - `autopvs1_db_str`: Name of dir created when `autopvs1_db` tar ball in decompressed. Default: `data`
 - `intervar_db_str`: Name of dir created when `intervar_db_str` tar ball in decompressed. Default: `intervardb`

### InterVar Classification Workflow
This workflow is a critical component in generating scoring metrics needed to classify pathogenicity of variants.
The workflow has three steps:
1. Convert VCF input into ANNOVAR txt format
1. Run ANNOVAR to gather needed scoring metrics
1. Run InterVar

### AutoPVS1
An additional pathogenicty scoring tool.
Documentation for this can be found [here](autopvs1/README.md)