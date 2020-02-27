FROM ubuntu:latest

MAINTAINER Kenny Westerman <kewesterman@mgh.harvard.edu>

RUN apt-get update && apt-get install -y wget unzip


RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip \
	&& unzip plink2_linux_avx2.zip


#COPY format_probabel_phenos.py format_probabel_output.py /

