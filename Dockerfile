FROM ubuntu

MAINTAINER Hendrik Feuerstein (hfeuerstein)

#Install without question 
ENV DEBIAN_FRONTEND noninteractive

# ARG and ENV for installing miniconda3
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update 
RUN apt-get install nano

# Install wget command 
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

#Install miniconda3 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py39_4.11.0-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py39_4.11.0-Linux-x86_64.sh 

# Set conda channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda 
RUN conda config --add channels conda-forge

# Installation of packages for bash script: 
RUN conda install fastqc

RUN conda install trimmomatic 

RUN conda install star

RUN conda install python=3.7.8

RUN conda install htseq

RUN conda install samtools=1.9


