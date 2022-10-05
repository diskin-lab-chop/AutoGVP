FROM ubuntu:20.04
LABEL maintainer="Miguel Browmn <brownm28@chop.edu>"
LABEL description="Docker file with modified version of autopvs1"

ENV D3B_AUTOPVS1_VERSION=v0.2.0

RUN apt update -y && apt install -y wget python3 pigz && \
wget https://github.com/diskin-lab-chop/pathogenicity-assessment/archive/refs/tags/${D3B_AUTOPVS1_VERSION}.tar.gz && \
rm -rf data cwl && \
pip3 install pysam==0.17.0 pyfaidx && \
apt remove -y wget