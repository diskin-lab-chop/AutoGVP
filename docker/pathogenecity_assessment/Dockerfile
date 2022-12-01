FROM rocker/tidyverse:4.2
MAINTAINER naqvia@chop.edu
WORKDIR /rocker-build/
RUN RSPM="https://packagemanager.rstudio.com/cran/2022-10-07" \
  && echo "options(repos = c(CRAN='$RSPM'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

COPY scripts/install_bioc.r .
COPY scripts/install_github.r .

# install wget, python and perl
RUN apt update && apt install -y perl pigz python3 wget

# Install apt-getable packages to start
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Install dev libraries and curl
RUN apt update && apt install -y \
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

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
  default-jdk \
  libxt6

# install R packages
RUN ./install_bioc.r \
    Biobase \
    BiocManager \
    optparse \
    vroom

ADD Dockerfile .
