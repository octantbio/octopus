# Start with rocker/tidyverse base image
FROM ubuntu:18.04

#===============================================================================
# Enviornment Variables

ENV bbmap_version 38.41
ENV PATH /usr/local/bbmap:${PATH}

ENV fgbio_version 1.0.0

ENV mlr_version 5.4.0

ENV minimap_version 2.16

ENV htslib_version 1.9
ENV bcftools_version 1.9
ENV samtools_version 1.9

ENV spades_version 3.13.0
ENV PATH /usr/local/SPAdes-${spades_version}-Linux/bin:${PATH}

ENV freebayes_version 1.3.1

ENV starcode_version 1.3

#===============================================================================
# Install extra *nix utils

RUN apt-get update \
    && apt-get install -y \
    pigz \
    vim \
    less \
    curl \
    wget \
    parallel \
    openjdk-8-jdk \
    python3-dev

RUN echo 'will cite' | parallel --citation 1> /dev/null 2> /dev/null
RUN ln -s /usr/bin/python3 /usr/bin/python

#===============================================================================
# Local Programs

# BBMap
WORKDIR /usr/local
RUN set -u; \
    wget -O bbmap.tar.gz "https://sourceforge.net/projects/bbmap/files/BBMap_${bbmap_version}.tar.gz/download" \
    && tar xzf bbmap.tar.gz \
    && rm bbmap.tar.gz

# FGBio
WORKDIR /usr/local/bin
RUN set -u; \
    wget "https://github.com/fulcrumgenomics/fgbio/releases/download/${fgbio_version}/fgbio-${fgbio_version}.jar" \
    && chmod 0644 /usr/local/bin/fgbio-${fgbio_version}.jar \
    && mv /usr/local/bin/fgbio-${fgbio_version}.jar /usr/local/bin/fgbio.jar \
    && echo '#!/bin/bash\njava -jar /usr/local/bin/fgbio.jar "$@"' > /usr/local/bin/fgbio \
    && chmod a+x /usr/local/bin/fgbio

# mlr
WORKDIR /usr/local/bin
RUN set -u; \
    wget "https://github.com/johnkerl/miller/releases/download/${mlr_version}/mlr.linux.x86_64" \
    && chmod 0777 mlr.linux.x86_64 \
    && mv mlr.linux.x86_64 mlr

# minimap2
WORKDIR /usr/local
RUN set -u; \
    wget -O minimap2.tar.bz2 "https://github.com/lh3/minimap2/releases/download/v${minimap_version}/minimap2-${minimap_version}_x64-linux.tar.bz2" \
    && tar xjf minimap2.tar.bz2 \
    && ln -s "/usr/local/minimap2-${minimap_version}_x64-linux/minimap2" /usr/local/bin/minimap2 \
    && rm minimap2.tar.bz2

# htslib
WORKDIR /usr/local
RUN apt-get install -y \
    build-essential \
    automake \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev

RUN set -u; \
    wget -O htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2" \
    && tar xjf htslib.tar.bz2 \
    && cd "htslib-${htslib_version}" \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm htslib.tar.bz2

# bcftools
WORKDIR /usr/local
RUN set -u; \
    wget -O bcftools.tar.bz2 "https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2" \
    && tar xjf bcftools.tar.bz2 \
    && cd "bcftools-${bcftools_version}" \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm bcftools.tar.bz2

# samtools
WORKDIR /usr/local
RUN set -u; \
    wget -O samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2" \
    && tar xjf samtools.tar.bz2 \
    && cd "samtools-${samtools_version}" \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm samtools.tar.bz2

# SPADES
WORKDIR /usr/local
RUN set -u; \
    wget "https://github.com/ablab/spades/releases/download/v${spades_version}/SPAdes-${spades_version}-Linux.tar.gz" \
    && tar xzf "SPAdes-${spades_version}-Linux.tar.gz" \
    && rm "SPAdes-${spades_version}-Linux.tar.gz" \
    && apt-get install -y python3-distutils

# freebayes
WORKDIR /usr/local/bin
RUN set -u; \
    wget "https://github.com/ekg/freebayes/releases/download/v${freebayes_version}/freebayes-v${freebayes_version}" \
    && mv "freebayes-v${freebayes_version}" freebayes \
    && chmod a+x freebayes

# starcode
WORKDIR /usr/local
RUN set -u; \
    wget -O starcode.tar.gz "https://github.com/gui11aume/starcode/archive/${starcode_version}.tar.gz" \
    && tar xzf starcode.tar.gz \
    && cd "starcode-${starcode_version}" \
    && make \
    && cp starcode /usr/local/bin/starcode \
    && rm ../starcode.tar.gz

# python requirements 
WORKDIR /usr/local
RUN set -u; \
    curl https://bootstrap.pypa.io/pip/3.6/get-pip.py | python3 \
    && pip install \
    pandas \
    numpy \
    stringdist \
    requests \
    pytest \
    openpyxl==3.0.7

#===============================================================================
# R install - lifted from
# https://github.com/rocker-org/rocker/blob/master/r-apt/bionic/Dockerfile
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        littler \
        r-base \
        r-base-dev \
        r-recommended \
        r-cran-rcpp \
    && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

# tidyverse install. adopted from
# https://github.com/rocker-org/rocker-versioned2/blob/master/scripts/install_tidyverse.sh
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    unixodbc-dev \
    && rm -rf /var/lib/apt/lists/*

# defaulted ncpus to 16
RUN install2.r --error \
    tidyverse \
    devtools \
    rmarkdown \
    vroom \
    gert \
    optparse


#===============================================================================
# set up this version of the pipeline inside the image under /opt
RUN mkdir -p /opt/octopus/data /opt/octopus/pipeline
# copy latest source code, test data, and Makefile
COPY test/ /opt/octopus/test
COPY src/ /opt/octopus/src/
COPY LICENSE Makefile /opt/octopus/

# Really important: let Python print emojis
# (with LANG unset python defaults to ASCII and we get UnicodeEncodeError)
ENV PYTHONIOENCODING utf8

# Create a directory that allows anyone to dump anything into it.
# It doesn't matter where this directory is, it just has to exist.
RUN mkdir /data && chmod 777 /data
WORKDIR /data

ENTRYPOINT ["/opt/octopus/src/main.py"]
