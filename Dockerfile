FROM rocker/tidyverse:4.2
MAINTAINER naqvia@chop.edu
WORKDIR /rocker-build/
RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

ENV D3B_AUTOPVS1_VERSION=v0.2.0a

COPY scripts/install_bioc.r .
COPY scripts/install_github.r .
COPY ./autopvs1 /autopvs1

# install wget, python and perl
RUN apt update -y && apt install -y perl pigz python3 wget python3-pip libbz2-dev liblzma-dev libxt-dev libproj-dev libv8-dev cpanminus libgdal-dev libgmp-dev libmpfr-dev

# Install apt-getable packages to start
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Install dev libraries and curl
RUN apt update -y && apt install -y \
    build-essential \
    bzip2 \
    cpanminus \
    cmake \
    curl \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgdal-dev \
    libgmp-dev \
    liblzma-dev \
    libmpfr-dev \
    libncurses5-dev \
    libproj-dev \
    libreadline-dev \
    libssl-dev \
    libv8-dev \
    libxt-dev \
    zlib1g-dev

## install python modules
RUN pip3 install pysam pyfaidx

# install R packages
RUN ./install_bioc.r \
    Biobase \
    BiocManager \
    optparse \
    vroom

# specific to autopvs1
#RUN rm -rf data cwl

ADD Dockerfile .
