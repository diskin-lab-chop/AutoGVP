FROM ubuntu:20.04
LABEL maintainer="Miguel Browmn <brownm28@chop.edu>"
LABEL description="Docker file with modified version of autopvs1"

ENV D3B_AUTOPVS1_VERSION=v0.2.0a

RUN apt update -y && apt install -y wget python3 python3-pip pigz
COPY ./autopvs1 /autopvs1
RUN rm -rf data cwl
RUN pip3 install pysam==0.17.0 pyfaidx && \
apt remove -y wget