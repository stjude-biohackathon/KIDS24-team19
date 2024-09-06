#! Dockerfile for installing tools and copying files and scripts for main CoCoRV package

# Pull base image.
FROM ubuntu:22.04

# Set maintainer.
LABEL  description='Dockerfile for installing python 3.10, tools and copying files and scripts for main CoCoRV package' \
    maintainer="Saima Sultana Tithi, St. Jude Children's Research Hospital, <saimasultana.tithi@stjude.org>"

# Define working directory.
WORKDIR /opt

RUN cd /opt

RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common \
    dirmngr libssl-dev libcurl4-openssl-dev libxml2-dev libopenblas-dev libopenblas-base liblapack3 wget curl bzip2 build-essential zlib1g-dev libbz2-dev liblzma-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -vxjf htslib-1.20.tar.bz2 && cd htslib-1.20 && ./configure && make && make install

RUN apt-get update && apt-get install -y \
    tabix bcftools bedtools \
    openjdk-8-jre-headless \
    g++ \
    python3.10 python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get install -y python3.10-dev libpq-dev && \
    python3.10 -m pip install --trusted-host pypi.python.org pip --upgrade && \
    python3.10 -m pip install --trusted-host pypi.python.org psycopg2-binary Cython numpy scikit-learn gnomad==0.7.1 hail==0.2.128 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PYSPARK_SUBMIT_ARGS="--driver-memory 40G pyspark-shell"

# Copy Dockerfile and Changelog.
COPY Dockerfile /opt/
#COPY Changelog /opt/

COPY utilities /opt/utilities
RUN chmod -R a+x /opt/utilities

# by default /bin/bash is executed
CMD ["/bin/bash"]
