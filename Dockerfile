FROM rocker/tidyverse:4.4.0
LABEL maintainer="Ryan Corbett (rcorbett@childrensnational.org)"
WORKDIR /rocker-build/

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl 
    
## install wget
RUN apt update -y && apt install -y wget bzip2 libbz2-dev

RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2  && \
    tar -xvjf bcftools-1.17.tar.bz2 && rm -f bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    make && mv /rocker-build/bcftools-1.17/bcftools /bin/.

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19', ask = FALSE)"

# install R packages
RUN R -e 'BiocManager::install(c( \
    "Biobase", \
    "lubridate", \
    "optparse", \
    "vcfR" \
))' 

RUN R -e "install.packages('pak', dependencies=TRUE); pak::pkg_install('url::https://github.com/tidyverse/vroom/archive/refs/tags/v1.7.1.tar.gz', ask = FALSE); remove.packages('pak')"

# AutoGVP
COPY . /rocker-build/AutoGVP/

ADD Dockerfile .
