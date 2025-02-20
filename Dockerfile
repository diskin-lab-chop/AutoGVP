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
    "BiocManager", \
    "lubridate", \
    "optparse", \
    "vcfR", \
    "vroom" \
))' 
    
# AutoGVP
RUN git clone https://github.com/diskin-lab-chop/AutoGVP.git

ADD Dockerfile .
