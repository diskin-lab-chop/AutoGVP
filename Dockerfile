FROM rocker/tidyverse:4.2
MAINTAINER naqvia@chop.edu
WORKDIR /rocker-build/
RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

COPY scripts/install_bioc.r .
COPY scripts/install_github.r .

## install wget
-RUN apt update -y && apt install -y wget bzip2

# install R packages
RUN ./install_bioc.r \
    Biobase \
    BiocManager \
    optparse \
    vroom

ADD Dockerfile .
