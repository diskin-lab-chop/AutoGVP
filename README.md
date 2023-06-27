# AutoGVP: Automated Germline Variant Pathogenicity
Jung Kim, Ammar S. Naqvi, Rebecca Kaufman, Zalman Vaksman, Miguel A. Brown, Ryan J. Corbett, Daniel P. Miller, Zhuangzhuang Geng, Phillip B. Storm, Adam C. Resnick, Jo Lynne Rokita, Douglas R. Stewart, Sharon J. Diskin

## AutoGVP Workflow  
<img src="https://github.com/diskin-lab-chop/pathogenicity-assessment/blob/b461f6248ea3bd472d646d3dd39445b616fa9295/figures/germline-pathogenecity_flow.png" align="center" width = "600">

## Clone the AutoGVP repository
```bash
git clone https://github.com/diskin-lab-chop/AutoGVP.git
```

## Docker set-up

### docker pull
```bash
docker pull pgc-images.sbgenomics.com/naqvia/autogvp:latest
```
cd to your clone of `AutoGVP`

### docker run
Replace <CONTAINER_NAME> with any name and run the command below:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -v $PWD:/home/rstudio/AutoGVP pgc-images.sbgenomics.com/naqvia/autogvp:latest
```
### docker execute
```bash
docker exec -ti <CONTAINER_NAME> bash
```

### cd to AutoGVP directory within docker
```bash
cd /home/rstudio/AutoGVP
```

## How to Run

### CAVATICA input ###
1. Run the [Kids First Germline Annotation Workflow](https://github.com/kids-first/kf-germline-workflow/blob/v0.4.4/docs/GERMLINE_SNV_ANNOT_README.md) first.
This workflow currently annotates variants with ClinVar (2022-05-07).
2. Run the [Pathogenicity Preprocessing Workflow](https://github.com/d3b-center/D3b-Pathogenicity-Preprocessing), which performs ANNOVAR with InterVar and AutoPVS1 annotations.
3. Run AutoGVP.

AutoGVP Requirements (recommended to place all in the `input/` folder):
- VEP-, ANNOVAR-, and ClinVar annotated VCF file (`*VEP.vcf`)
- ANNOVAR multianno file (`*hg38_multianno.txt`)
- InterVar file (`*intervar.hg38_multianno.txt.intervar`)
- AutoPVS1 file (`*autopvs1.txt`)
- `variant_summary.txt`
- `submission_summary.txt` (can retrieve with `download_db_files.sh`)
Note: the variant_summary and submission_summary files need to be uncompressed.

Example run:

```bash
bash download_db_files.sh
```
```r
Rscript 01-annotate_variants_custom_input.R --vcf <*.vcf> --multianno <*multianno.txt> --intervar <*hg38_multianno.txt.intervar> --autopvs1 <*autopvs1.txt --output <output_prefex> --submission input/variant_summary.txt --submission_summary input/submission_summary.txt
```

### Custom (non-CAVATICA) input ###
1. Annotate the germline VCF with [VEP](https://github.com/Ensembl/ensembl-vep).
Note: It is recommended to run VEP 104 to ensure optimal tool compatibility since AutoPVS1 hg38 uses gene symbols from VEP 104.
Alternatively, if using VEP > 104, it is recommended to lift over the gene symbols in the `PVS1.level` file located in the AutoPVS1 data folder using this [custom python script](https://github.com/d3b-center/D3b-DGD-Collaboration/blob/main/scripts/update_gene_symbols.py) where `hgnc_tsv` is the gene name database TSV file from the monthly HGNC server [here](https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/).
Example command, with results used to replace `PVS1.level` file. 
```python
python3 D3b-DGD-Collaboration/scripts/update_gene_symbols.py -g hgnc_complete_set_2021-06-01.txt -f PVS1.level -z GENE level -u GENE -o results --explode_records 2> old_new.log
```

2. Run [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) with the following options (to create the VCF input for AutoGVP):
```perl
perl table_annovar.pl input/test_hg38_selected_VEP_annotated.vcf hg38 --buildver hg38 --out test_hg38_selected --remove --protocol gnomad211_exome,gnomad211_genome --operation f,f --vcfinput
```
3. Run [InterVar](https://github.com/WGLab/InterVar) with the following command:
```python
python InterVar.py -b hg38 -i input/test_hg38_selected_VEP_annotated.vcf --input_type=VCF -o test_hg38_selected
```
4. Run [AutoPVS1 v2.0](https://github.com/JiguangPeng/autopvs1/releases/tag/v2.0).
5.  Optional: provide a ClinVar VCF file. If not supplied by the user, the most recent ClinVar file will be downloaded with `download_db_files.sh` and used in AutoGVP.
6. Run AutoGVP.

AutoGVP Requirements (recommended to place all in the `input/` folder):
- VEP-annotated VCF (`*VEP.vcf`)
- ANNOVAR file (`*hg38_multianno.txt`)
- InterVar file (`*intervar.hg38_multianno.txt.intervar`)
- AutoPVS1 file (`*autopvs1.txt`)
- ClinVar VCF (`clinvar_yyyymmdd.vcf.gz` optional user input or `clinvar.vcf.gz` will be downloaded with `download_db_files.sh`)
- `variant_summary.txt`
- `submission_summary.txt` (can retrieve with `download_db_files.sh`)
Note: the variant_summary and submission_summary files need to be uncompressed.

Example run:
```bash
bash download_db_files.sh
```
```r
Rscript 01-annotate_variants_custom_input.R --vcf input/testing_010423_VEP.vcf --multianno input/testing_010423.hg38_multianno.txt --intervar input/testing_010423_intervar.hg38_multianno.txt.intervar --autopvs1 input/testing_010423_autopvs1.txt --clinvar input/clinvar_20211225.vcf.gz --output SRRT0182 --submission input/variant_summary.txt --submission_summary input/submission_summary.txt
```

### Step by step workflow
1. Read in ClinVar file
2. Create `vcf_id` (CHR, POS, REF, ALT)
3. Annotate stars based on `CLNREVSTAT`*
4. Report and save significant call if Stars are not 0 or 1 and not B/P/LB/LP/VUS
5. Identify variants that are ambiguous (`criteria_provided,_conflicting_interpretations`), generate table to then check against submission file
6. Identify variants that need further annotations or possible re-adjustment (Stars 0 or 1)
7. Retrieve and store interVar results file into table and create vcf_id
8. Retrieve and store corresponding autopvs1 results file into table and create vcf_id
9. Merge interVar and autopvs1 tables by matching vcf_ids
10. Create columns for `evidencePVS1`, `evidencePS`, `evidencePM`, `evidencePP`, `evidenceBP`, `evidencePM` and `evidenceBA1` (variables that may need re-adjusting) by parsing InterVar: InterVar and Evidence column
11. Indicate if there needs to be recalculation (if `evidencePVS1` == 1)
12. If not, take interVar significant call as final_call
13. Go through entries and adjust evidence variables above based on criterion*
14. Report final call based on recalculated evidence variables

#### Annotating Stars
```
1 = 'criteria_provided,_single_submitter','Needs resolving'
2 = 'criteria_provided,_multiple_submitters'
3 = 'reviewed_by_expert_panel'
4 = 'practice_guideline'
0 = 'criteria_provided,_conflicting_interpretations'
```

#### Re-calculation adjustments
```
if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
if criterion is NF6|SS6|DEL7|DEL11|IC3 then PVS1 = 0; PP = PP+1;
if criterion is IC4 then PVS1 = 0; PP = PP+1;
if criterion is na|NF0  then PVS1 = 0;

New ClinSig
Pathogenic - Criteria 1
  (i) 1 Very strong (PVS1) AND
        (a) ≥1 Strong (PS1–PS4) OR
        (b) ≥2 Moderate (PM1–PM6) OR
        (c) 1 Moderate (PM1–PM6) and 1 supporting (PP1–PP5) OR
        (d) ≥2 Supporting (PP1–PP5)

Pathogenic - Criteria 2
  (ii) ≥2 Strong (PS1–PS4) OR

Pathogenic - Criteria 3
  (iii) 1 Strong (PS1–PS4) AND
        (a)≥3 Moderate (PM1–PM6) OR
        (b)2 Moderate (PM1–PM6) AND ≥2 Supporting (PP1–PP5) OR
        (c)1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5)

Likely pathogenic
        (i) 1 Very strong (PVS1) AND 1 moderate (PM1– PM6) OR
        (ii) 1 Strong (PS1–PS4) AND 1–2 moderate (PM1–PM6) OR
        (iii) 1 Strong (PS1–PS4) AND ≥2 supporting (PP1–PP5) OR
        (iv)  ≥3 Moderate (PM1–PM6) OR
        (v) 2 Moderate (PM1–PM6) AND ≥2 supporting (PP1–PP5) OR
        (vi) 1 Moderate (PM1–PM6) AND ≥4 supporting (PP1–PP5)

Benign
        (i) 1 Stand-alone (BA1) OR
        (ii) ≥2 Strong (BS1–BS4)

Likely Benign
        (i) 1 Strong (BS1–BS4) and 1 supporting (BP1– BP7) OR
        (ii) ≥2 Supporting (BP1–BP7)

Uncertain  significance
        (i) non of the criteria were met.
        (ii) Benign and pathogenic are contradictory.
```

## Code Author

Ammar S. Naqvi ([@naqvia](https://github.com/naqvia))

## Contact

For questions, please submit an issue or send an email to Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)): rokita@chop.edu
