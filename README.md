# AutoGVP: Automated Germline Variant Pathogenicity
Jung Kim, Ammar S. Naqvi, Ryan J. Corbett, Rebecca Kaufman, Zalman Vaksman, Miguel A. Brown,  Daniel P. Miller, Zhuangzhuang Geng, Phillip B. Storm, Adam C. Resnick, Jo Lynne Rokita, Douglas R. Stewart, Sharon J. Diskin

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
3. Download required db files, including variant submissions files (`variant_summary.txt` and `submission_summary.txt`) using download_db_files.sh (see below).
4. Run `select-clinVar-submissions.R` within `AutoGVP/` to consolidate ClinVar submission calls prior to running AutoGVP wrapper script. This step can be ran once and the `ClinVar-selected-submissions.tsv` file could be used for subsequent runs if clinVar versions remain the same. 
5. Run AutoGVP wrapper script using `--workflow="cavatica"`

AutoGVP Requirements (recommended to place all in the `input/` folder):
- VEP-, ANNOVAR-, and ClinVar annotated VCF file (`*VEP.vcf`)
- ANNOVAR multianno file (`*hg38_multianno.txt`)
- InterVar file (`*intervar.hg38_multianno.txt.intervar`)
- AutoPVS1 file (`*autopvs1.txt`)
- Variant submissions file (`ClinVar-selected-submissions.tsv` generated by `select-clinVar-submissions.R`)


Example run:

Download ClinVar files:

```bash
bash download_db_files.sh
```

Run `select-clinVar-submissions.R`:

```
Rscript `select-clinVar-submissions.R` --variant_summary <variant_summary.txt.gz> --submission_summary <submission_summary.txt.gz>
```

Run AutoGVP:

```r
bash run_autogvp.sh --workflow="cavatica" \
--vcf=input/test_pbta.single.vqsr.filtered.vep_105.vcf \
--filter_criteria=<filter criteria> \
--intervar=input/test_pbta.hg38_multianno.txt.intervar \
--multianno=input/test_pbta.hg38_multianno.txt \
--autopvs1=input/test_pbta.autopvs1.tsv \
--outdir=../results \
--out="test_pbta"

```

### Custom (non-CAVATICA) input ###
1. Annotate the germline VCF with [VEP](https://github.com/Ensembl/ensembl-vep).
Note: It is recommended to run VEP 104 to ensure optimal tool compatibility since AutoPVS1 hg38 uses gene symbols from VEP 104.
Alternatively, if using VEP > 104, it is recommended to lift over the gene symbols in the `PVS1.level` file located in the AutoPVS1 data folder using this [custom python script](https://github.com/d3b-center/D3b-DGD-Collaboration/blob/main/scripts/update_gene_symbols.py) where `hgnc_tsv` is the gene name database TSV file from the monthly HGNC server [here](https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/).
Example command, with results used to replace `PVS1.level` file.
```python
python3 D3b-DGD-Collaboration/scripts/update_gene_symbols.py -g hgnc_complete_set_2021-06-01.txt -f PVS1.level -z GENE level -u GENE -o results --explode_records 2> old_new.log
```

2. Run [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) with the following options in order to create ANNOVAR annotated file using VCF input:
```perl
perl table_annovar.pl input/test_VEP.vcf hg38 --buildver hg38 --out test_VEP --remove --protocol gnomad211_exome,gnomad211_genome --operation f,f --vcfinput
```
3. Run [InterVar](https://github.com/WGLab/InterVar) with the following command:
```python
python InterVar.py -b hg38 -i input/test_VEP.vcf --input_type=VCF -o test_VEP
```
4. Run [AutoPVS1 v2.0](https://github.com/JiguangPeng/autopvs1/releases/tag/v2.0).
5. Optional: provide a ClinVar VCF file. If not supplied by the user, the most recent ClinVar file will be downloaded with `download_db_files.sh` and used in AutoGVP.
6. Download required db files, including variant submissions (`variant_summary.txt` and `submission_summary.txt`) and clinVar files using download_db_files.sh (see below).
7. Run `select-clinVar-submissions.R` within `AutoGVP/` to consolidate ClinVar submission calls prior to running AutoGVP wrapper script. This step can be ran once and the `ClinVar-selected-submissions.tsv` file could be used for subsequent runs if clinVar versions remain the same.
8. Run AutoGVP wrapper script using `--workflow="custom"`

AutoGVP Requirements (recommended to place all in the `input/` folder):
- VEP-annotated VCF (`*VEP.vcf`)
- ANNOVAR file (`*hg38_multianno.txt`)
- InterVar file (`*hg38_multianno.txt.intervar`)
- AutoPVS1 file (`*autopvs1.txt`)
- ClinVar VCF (`clinvar_yyyymmdd.vcf.gz` optional user input or `clinvar.vcf.gz` will be downloaded with `download_db_files.sh`)
- Variant submissions file (`ClinVar-selected-submissions.tsv` generated by `select-clinVar-submissions.R`)

Example run:

Download ClinVar files:

```bash
bash download_db_files.sh
```

Run `select-clinVar-submissions.R`:

```
Rscript `select-clinVar-submissions.R` --variant_summary <variant_summary.txt.gz> --submission_summary <submission_summary.txt.gz>
```

Run AutoGVP:
```r
bash run_autogvp.sh --workflow="custom" \
--vcf=input/test_VEP.vcf \
--filter_criteria=<filter criteria>
--clinvar=input/clinvar.vcf.gz \
--intervar=input/test_VEP.hg38_multianno.txt.intervar \
--multianno=input/test_VEP.vcf.hg38_multianno.txt \
--autopvs1=input/test_autopvs1.txt \
--outdir=../results \
--out="test_custom"
```

### Step by step workflow

* __Filter VCF file__. By default, AutoGVP filters based on `FILTER` column (`PASS` or `.`). Other criteria can be specified by `filter_criteria` argument as follows:

```
filter_criteria="INFO/AF>=0.2 INFO/DP>=15"
```

* __Run Pathogenicity Assessment__. The R scripts `01-annotate_variants_CAVATICA_input.R` and `01-annotate_variants_custom_input.R` perform the following steps:

  1. Read in ClinVar-annotated VCF file
  2. Assign ClinVar stars based on `CLNREVSTAT`*
  3. For ClinVar variants, report `CLINSIG` as final call; resolve ambiguous variants (`criteria_provided,_conflicting_interpretations`) by checking against ClinVar variant submission file
  4. Identify variants that need further Intervar annotation and possible re-adjustment (variants with 0 stars or not in ClinVar database)
  5. Load and merge ANNOVAR multianno, InterVar, and AutoPVS1 files
  6. Create columns for `evidencePVS1`, `evidencePS`, `evidencePM`, `evidencePP`, `evidenceBP`, `evidencePM` and `evidenceBA1` (variables that may need re-adjusting) by parsing `InterVar: InterVar and Evidence` column
  7. Adjust evidence columns based on AutoPVS1 `criterion` column
  8. Report InterVar final call (if unadjusted) or final call based on re-calculated evidence variables (if adjusted)
  9. Save output

* __Parse VCF file__. `parse_vcf.sh` converts the VCF file to a TSV file with INFO fields as tab-separated columns.

* __Resolve gene annotations and produce final output files__ (`04-filter_gene_annoations.R`)

  1. Read in parsed VCF file, and select single VEP annotation based on `PICK` column (`PICK == 1`).
  See R script for criteria used to select pick transcripts.
  2. Merge gene annotation with AutoGVP results
  3. Select columns to retain in final output, and save full and abridged output files


#### Annotating Stars (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/, https://www.nejm.org/doi/10.1056/NEJMsr1406261)
```
1 = 'criteria_provided,_single_submitter','criteria_provided,_conflicting_interpretations'
2 = 'criteria_provided,_multiple_submitters'
3 = 'reviewed_by_expert_panel'
4 = 'practice_guideline'
0 = 'no_assertion_provided','no_assertion_criteria_provided','no_assertion_for_the_individual_variant'
```

#### Re-calculation adjustments (https://onlinelibrary.wiley.com/doi/10.1002/humu.23626)
```
if criterion is NF1|SS1|DEL1|DEL2|DUP1|IC1 then PVS1=1
if criterion is NF3|NF5|SS3|SS5|SS8|SS10|DEL4|DEL8|DEL6|DEL10|DUP3|IC2 then PVS1 = 0; PS = PS+1
if criterion is NF6|SS6|SS9|DEL7|DEL11|IC3 then PVS1 = 0; PM = PM+1;
if criterion is IC4 then PVS1 = 0; PP = PP+1;
if criterion is na|NF0|NF2|NF4|SS2|SS4|SS7|DEL3|DEL5|DEL9|DUP2|DUP4|DUP5|IC5  then PVS1 = 0;

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

### Output

AutoGVP produces an abridged output file with minimal information needed to interpret variant pathogenicity, as well as a full output with >100 variant annotation columns.

#### Abridged output example:

chr | start | ref | alt | rs_id | gene_symbol_vep | variant_classification_vep | HGVSg | HGVSc | HGVSp | autogvp_call | autogvp_call_reason | clinvar_stars | clinvar_clinsig | intervar_evidence
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
chr1 | 1332490 | C | T | rs201607183 | TAS1R3 | missense_variant | chr1:g.1332490C>T | c.959C>T | p.Thr320Met | Uncertain_significance | ClinVar | 1 | Uncertain_significance | InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 0, 0, 0, 0, 0, 0] PP=[0, 0, 1, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]
chr1 | 1390349 | C | T | rs769726291 | CCNL2 | missense_variant | chr1:g.1390349C>T | c.887G>A | p.Gly296Asp | Uncertain_significance | InterVar | NA | NA | InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]

*NOTE: gnomAD v.3.1.1 non-cancer AF popmax values (`gnomad_3_1_1_AF_non_cancer`) will also be included in abridged output when provided.

See [here](https://github.com/diskin-lab-chop/AutoGVP/blob/main/AutoGVP/input/output_colnames.tsv) for list of columns included in full output


## Code Authors

Ammar S. Naqvi ([@naqvia](https://github.com/naqvia)) and Ryan J. Corbett ([@rjcorb](https://github.com/rjcorb))

## Contact

For questions, please submit an issue or send an email to Jo Lynne Rokita ([@jharenza](https://github.com/jharenza)): rokita@chop.edu
