FROM ubuntu:20.04
LABEL maintainer="Miguel Browmn <brownm28@chop.edu>"
LABEL description="Docker file with software necessary to run intervar, a pathogenicity predictor, and annovar, output of wihch is necessary for prediction"

ENV INTERVAR_VERSION=v2.2.1

RUN apt update && apt install -y wget python3 perl pigz \
&& wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz \
&& tar -xzf annovar.latest.tar.gz && rm annovar.latest.tar.gz \
&& wget https://raw.githubusercontent.com/WGLab/InterVar/${INTERVAR_VERSION}/Intervar.py \
&& wget https://raw.githubusercontent.com/WGLab/InterVar/${INTERVAR_VERSION}/config.ini \
&& apt autoremove -y