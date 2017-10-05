FROM ubuntu:xenial
MAINTAINER Christopher J. Yoon <cjyoon@wustl.edu>
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository ppa:jonathonf/python-3.6 && \
  apt-get update && apt-get install -y python3.6
RUN sed -i 's/archive.ubuntu.com/mirror.us.leaseweb.net/' /etc/apt/sources.list \
    && sed -i 's/deb-src/#deb-src/' /etc/apt/sources.list \
    && apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
    build-essential \
    ca-certificates \
    gcc \
    git \
    libpq-dev \
    make \
    pkg-config \
    python3-pip \
    aria2 \
    python3-tk \
    bwa \
    libboost-all-dev \
    libbz2-dev \
    liblzma-dev \
    samtools \
    bedtools \
    r-base \
    wget \
    && apt-get autoremove -y \
    && apt-get clean
    
# install java8
RUN apt-get update && \
  apt-get install -y software-properties-common python-software-properties
RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer
# download GATK jar
RUN mkdir /home/software && cd /home/software && wget https://github.com/broadinstitute/gatk/releases/download/4.beta.5/gatk-4.beta.5.zip && \
  unzip gatk-4.beta.5.zip
  
# Install R and necessary packages
RUN  cd /home/software/ &&  \
  wget https://bioconductor.org/packages/release/bioc/src/contrib/DNAcopy_1.50.1.tar.gz && \
  R CMD INSTALL DNAcopy_1.50.1.tar.gz  && \
  wget https://cran.r-project.org/src/contrib/naturalsort_0.1.3.tar.gz && \
  R CMD INSTALL naturalsort_0.1.3.tar.gz  && \
  wget https://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz && \
  R CMD INSTALL getopt_1.20.0.tar.gz  && \
  wget https://cran.rstudio.com/src/contrib/optparse_1.4.4.tar.gz && \
  R CMD INSTALL optparse_1.4.4.tar.gz 

# Download github repository for gatk4 CNV WXS pipeline
RUN cd /home/software/ && \
  git clone https://github.com/ding-lab/gatk4wxscnv.git
