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

## Get data files

<br>**Run shell script to get dummy test files**
```
bash download_data.sh
```
