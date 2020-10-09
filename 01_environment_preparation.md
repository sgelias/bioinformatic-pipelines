# Environment preparation

Aim: prepare a local or docker machine for bioinformatic analysis.

## Linux environment preparation (based in Linux Ubuntu >= 18.04 with root access)

## Install basics

```bash
apt update -y && apt upgrade -y
```

```bash
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
```

### Install the main dependencies for bioinformatics, as [cutadapt](https://pypi.org/project/cutadapt/), [biopython](https://pypi.org/project/biopython/), [ncbi-genome-download](https://pypi.org/project/ncbi-genome-download/), and [ncbi-acc-download](https://pypi.org/project/ncbi-acc-download/).

```bash
pip3 install --upgrade \
    cutadapt \
    biopython \
    ncbi-genome-download \
    ncbi-acc-download
```

### Get CPANminus from [pre-compiled source code](http://cpanmin.us/)

```bash
curl -L http://cpanmin.us | perl - --force Time::HiRes
```

### SRA Toolkit >= 2.9.6-1

```bash
wget \
    https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    && tar -xzf sratoolkit.current-ubuntu64.tar.gz \
    && rm sratoolkit.current-ubuntu64.tar.gz
```

```bash
ln -sfv /home/bioinfo/sratoolkit.2.9.6-1-ubuntu64/bin/* /usr/local/bin/
```

### Trimmomatic-0.39

```bash
wget \
    http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
    && unzip Trimmomatic-0.39.zip \
    && rm Trimmomatic-0.39.zip
```

### SPAdes-3.13

```bash
wget \
    http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz \
    && tar -xzf SPAdes-3.13.0-Linux.tar.gz \
    && rm SPAdes-3.13.0-Linux.tar.gz
```

### Quast-5.0.2

```bash
wget \
    https://sourceforge.net/projects/quast/files/quast-5.0.2.tar.gz \
    && tar -xzf quast-5.0.2.tar.gz \
    && rm quast-5.0.2.tar.gz
```

```bash
quast-5.0.2/install_full.sh && rm -rf quast_test_output
```

### Glimmer-3.02

```bash
wget \
    http://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz \
    && tar -xzf glimmer302b.tar.gz \
    && rm glimmer302b.tar.gz \
    && cd glimmer3.02/src/ \
    && make \
    && cd /home/bioinfo
```

