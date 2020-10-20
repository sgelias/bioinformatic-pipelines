FROM ubuntu:18.04

RUN apt update -y \
    && apt upgrade -y \
    && DEBIAN_FRONTEND="noninteractive" TZ="America" \
        apt install -y \
        wget \
        python3-pip \
        python3-dev \
        python-dev \
        zip \
        unzip \
        default-jdk \
        curl \
        nano \
        git \
        gcc \
        zlib1g-dev \
        zlib1g \
        build-essential \
        r-base \
        pkg-config \
        libfreetype6-dev \
        libpng-dev \
        python-matplotlib \
        ncbi-blast+ \
        bwa \
        hisat2 \
        bowtie2 \
        samtools

RUN pip3 install --upgrade \
    cutadapt \
    biopython \
    ncbi-genome-download \
    ncbi-acc-download

CMD [ "/bin/bash" ]
