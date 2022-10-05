# Pathogenicity assessment development

## InterVar Classification Workflow
This workflow is a critical component in generating scoring metrics needed to classify pathogenicity of variants.
It is recommended to have first run the [Kids First Germline Annotation Workflow](https://github.com/kids-first/kf-germline-workflow/blob/master/docs/GERMLINE_SNV_ANNOT_README.md) first.
The workflow has three steps:
1. Convert VCF input into ANNOVAR txt format
1. Run ANNOVAR to gather needed scoring metrics
1. Run InterVar

## AutoPVS1
An additional pathogenicty scoring tool.
Documentation for this can be found [here](autopvs1/README.md)