# germline-pathogenicity-assessment
Germline pathogenecity assessment

## Docker set-up

### docker pull and run
```
docker pull pgc-images.sbgenomics.com/naqvia/germline-pathogenicity-assessment:latest
docker run --platform linux/amd64 --name pathogenecity_anno -d -v $PWD:/home/rstudio/pathogenecity-assessment pgc-images.sbgenomics.com/naqvia/germline-pathogenicity-assessment:latest

```
### docker execute
```
docker exec -ti pathogenecity_anno bash
```

## How to Run
### dummy input data
```
input/*vcf
input/*intervar
input/*autopvs1.tsv
```

### Run annotation scripts (example)
```
bash run_annotator.sh -v input/BS_test_ad_hoc_genotyping.CGP.filtered.deNovo.vep.chr17_test.vcf -i input/BS_test_ad_hoc_annovar_humandb_hg38_intervar.hg38_multianno.chr17_test.txt.intervar -a input/BS_test_ad_hoc_annovar_humandb_hg38_intervar.chr17_test.autopvs1.tsv -w cavatica
```
