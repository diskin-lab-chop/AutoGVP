# AutoGVP: Automated Germline Variant Pathogenicity
[![DOI](https://zenodo.org/badge/525946788.svg)](https://zenodo.org/doi/10.5281/zenodo.10107711)

**This work has now been published**: [AutoGVP: a dockerized workflow integrating ClinVar and InterVar germline sequence variant classification.](https://doi.org/10.1093/bioinformatics/btae114)

Kim J^, Naqvi AS^, Corbett RJ, Kaufman RS, Vaksman Z, Brown MA, Miller DP, Phul S, Geng Z, Storm PB, Resnick AC, Stewart DR, Rokita JL+, Diskin SJ+. AutoGVP: a dockerized workflow integrating ClinVar and InterVar germline sequence variant classification. Bioinformatics. 2024 Mar 4;40(3):btae114. doi: 10.1093/bioinformatics/btae114. PMID: 38426335; PMCID: PMC10955249.

^Equal first authorship 
+Equal senior authorship

## AutoGVP Workflow  
<img src="https://github.com/diskin-lab-chop/AutoGVP/blob/main/figures/germline-pathogenecity_flow.png" align="center" width = "600">

For more detailed instructions, please visit the [user guide](https://github.com/diskin-lab-chop/AutoGVP/wiki/User-Guide) on our wiki.

## Clone the AutoGVP repository
```bash
git clone git@github.com:diskin-lab-chop/AutoGVP.git
```

## Docker set-up
1. Pull the docker image.
```bash
docker pull pgc-images.sbgenomics.com/diskin-lab/autogvp:v1.0.3
```
2. Navigate to the `AutoGVP` root directory
```bash
cd AutoGVP
```   
3. Start a docker image. Replace <CONTAINER_NAME> with any name and run the commands below:
```bash
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -v $PWD:/home/rstudio/AutoGVP pgc-images.sbgenomics.com/diskin-lab/autogvp:v1.0.3
docker exec -ti <CONTAINER_NAME> bash
```
4. Navigate to AutoGVP directory within the docker image
```bash
cd /home/rstudio/AutoGVP
```
5. Run AutoGVP (see example commands below).

## Dependencies
[VEP (v104)](https://github.com/Ensembl/ensembl-vep)<br>
[InterVar](https://github.com/WGLab/InterVar) <br>
[ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) <br>
[AutoPVS1 (v2.0)](https://github.com/JiguangPeng/autopvs1/releases/tag/v2.0) <br>
[bcftools (v1.17)](https://github.com/samtools/bcftools/) <br>

## How to Run AutoGVP
AutoGVP Requirements (recommended to place all in the `data/` folder):
- VEP-annotated VCF file (`*VEP.vcf`) or VEP- and ClinVar-annotated VCF file (CAVATICA workflow only). For CAVATICA workflow, AutoGVP will use ClinVar annotation from sample VCF when external ClinVar VCF file is not provided.   
- ANNOVAR multianno file (`*hg38_multianno.txt`)
- InterVar file (`*intervar.hg38_multianno.txt.intervar`)
- AutoPVS1 file (`*autopvs1.txt`)
- Variant submissions file (`ClinVar-selected-submissions.tsv` generated by `select-clinVar-submissions.R`)
- ClinVar VCF (`clinvar_yyyymmdd.vcf.gz` optional user input or `clinvar.vcf.gz` will be downloaded with `download_db_files.sh`). This is an optional input for CAVATICA workflow; if not provided, AutoGVP will expect ClinVar annotation in VEP-annotated sample VCF (see above).  

### Custom workflow example run
1. [Prepare input files](https://github.com/diskin-lab-chop/AutoGVP/wiki/User-Guide#custom-input-workflow---step-by-step) by running VEP, ANNOVAR, InterVar, and AutoPVS1.
2. Download [database files](https://github.com/diskin-lab-chop/AutoGVP/blob/main/scripts/download_db_files.sh): 
```
bash scripts/download_db_files.sh
```
3. Run `select-clinVar-submissions.R`. To customize conflicting interpretation resolution, users can provide a ClinGen Concept ID list to filter submissions against (`--conceptID_list`). When a list is provided, users can also determine how unsettled conflicts are resolved with the `--conflict_res` argument (`"latest"` or `"most_severe"`). 
For more details, see the [FAQ](https://github.com/diskin-lab-chop/AutoGVP/wiki/FAQ#how-can-i-create-my-own-concept-id-list).
Example command:
```
Rscript scripts/select-clinVar-submissions.R --variant_summary data/variant_summary.txt.gz --submission_summary data/submission_summary.txt.gz --outdir results --conceptID_list data/clinvar_cpg_concept_ids.txt --conflict_res "latest"
```
4. Run AutoGVP; if output of scripts/select-clinVar-submissions.R is not provided, the script will be run prior to starting pathogenicity assessment
```r
bash run_autogvp.sh --workflow="custom" \
--vcf=data/test_VEP.vcf \
--filter_criteria=<filter criteria>
--clinvar=data/clinvar.vcf.gz \
--intervar=data/test_VEP.hg38_multianno.txt.intervar \
--multianno=data/test_VEP.vcf.hg38_multianno.txt \
--autopvs1=data/test_autopvs1.txt \
--outdir=results \
--out="test_custom" \
--selected_clinvar_submissions=results/ClinVar-selected-submissions.tsv \
--variant_summary=data/variant_summary.txt.gz \
--submission_summary=data/submission_summary.txt.gz \
--conceptIDs=data/clinvar_cpg_concept_ids.txt \
--conflict_res="latest"
```

### CAVATICA workflow example run
1. Download [database files](https://github.com/diskin-lab-chop/AutoGVP/blob/main/scripts/download_db_files.sh): 
```
bash scripts/download_db_files.sh
```
2. Run `select-clinVar-submissions.R` (See custom workflow step 2 for optional conflict resolution parameters).
For more details, see the [FAQ](https://github.com/diskin-lab-chop/AutoGVP/wiki/FAQ#how-can-i-create-my-own-concept-id-list).
Example command:
```
Rscript scripts/select-clinVar-submissions.R --variant_summary data/variant_summary.txt.gz --submission_summary data/submission_summary.txt.gz --outdir results --conceptID_list data/clinvar_cpg_concept_ids.txt --conflict_res "latest"
```
2. Run AutoGVP; if output of scripts/select-clinVar-submissions.R is not provided, the script will be run prior to starting pathogenicity assessment

```r
bash run_autogvp.sh --workflow="cavatica" \
--vcf=data/test_pbta.single.vqsr.filtered.vep_105.vcf \
--filter_criteria=<filter criteria> \
--intervar=data/test_pbta.hg38_multianno.txt.intervar \
--multianno=data/test_pbta.hg38_multianno.txt \
--autopvs1=data/test_pbta.autopvs1.tsv \
--outdir=results \
--out="test_pbta" \
--selected_clinvar_submissions=results/ClinVar-selected-submissions.tsv \
--variant_summary=data/variant_summary.txt.gz \
--submission_summary=data/submission_summary.txt.gz \
--conceptIDs=data/clinvar_cpg_concept_ids.txt \
--conflict_res="latest"
```

### AutoGVP Output
AutoGVP produces an abridged output file with minimal information needed to interpret variant pathogenicity, as well as a full output with >100 variant annotation columns.

#### Abridged output example:

chr | start | ref | alt | rs_id | gene_symbol_vep | variant_classification_vep | HGVSg | HGVSc | HGVSp | autogvp_call | autogvp_call_reason | clinvar_stars | clinvar_clinsig | intervar_evidence
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
chr1 | 1332490 | C | T | rs201607183 | TAS1R3 | missense_variant | chr1:g.1332490C>T | c.959C>T | p.Thr320Met | Uncertain_significance | ClinVar | 1 | Uncertain_significance | InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 0, 0, 0, 0, 0, 0] PP=[0, 0, 1, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]
chr1 | 1390349 | C | T | rs769726291 | CCNL2 | missense_variant | chr1:g.1390349C>T | c.887G>A | p.Gly296Asp | Uncertain_significance | InterVar | NA | NA | InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0, 0] PM=[1, 1, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]

*NOTE: gnomAD v.3.1.1 non-cancer AF popmax values (`gnomad_3_1_1_AF_non_cancer`) will also be included in abridged output when provided.

See [here](https://github.com/diskin-lab-chop/AutoGVP/blob/main/data/output_colnames.tsv) for list of columns included in full output.


## Code Authors

Ammar S. Naqvi ([@naqvia](https://github.com/naqvia)) and Ryan J. Corbett ([@rjcorb](https://github.com/rjcorb))

## Contact

For questions, please submit an issue
