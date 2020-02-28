FROM ubuntu:latest

MAINTAINER Kenny Westerman <kewesterman@mgh.harvard.edu>

RUN apt-get update && apt-get install -y wget python3 python3-pip unzip
RUN pip3 install pandas scipy


RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip \
	&& unzip plink2_linux_avx2.zip


COPY format_plink2_phenos.py /

