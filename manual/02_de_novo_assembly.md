# De novo assembly (NOT FINISHED)

Aim: perform a genome assembly without reference sequence.

This tutorial explains how to:

1. Pull a docker environment (the same environment can be created using the steps described in **Environment preparation**) containing all tools necessary to perform the genome assembly and annotation;

2. Get the genome of *Enterobacter cloacae* subsp. *cloacae* (NCTC 9394 isolate) from Genbank;

3. Get the MEGAres database used for the prediction of resistance genes;

4. Cleaning the SRA file (filtering step) using Trimmomatic;

5. Sequence assembly using SPAdes;

6. Generate descriptive statistics for the assembly using Quast;

7. Predict resistance genes based on MEGAres database;

8. And finally, annotate.

### Pull docker image using pure docker or docker-compose (optional)

The docker image should be pulled from Docker Hub using the command below. Alternatively, a non-docker Linux machine can be prepared following the steps described in the section Environment preparation.

```bash
docker pull waldeyr/bioinfo:drug_resistance
```

Then, turn on the docker container using:

```bash
docker run --rm -it waldeyr/bioinfo:drug_resistance 
```

You can do it thought a docker-composer. Thus, create a .yml file:

```bash
touch docker-compose.yml && vim docker-compose.yml
```

and insert the content:

```yml
version: "3"

volumes:
    analysis: {}

services:
    machine:
        container_name: "drug_resistance"
        image: waldeyr/bioinfo:drug_resistance
        volumes:
            - analysis:/home/
```

Then, start the container using:

```bash
docker-compose run --rm machine
```

Using the docker-compose you have several advantages, including persist all parameters or/and create once all containers that composes the *orchestra*.

### Downloading Sequences

Running the container, create all directories needed in this tutorial.

```bash
mkdir /home/MEGAres /home/blast /home/data /home/glimmer
```

Then, download Illumina sequencing reads library from Genbank/SRA. The `--split-files` parameter split the paired-end file into forward and reverse files, respectively. The `fastq-dump` function are included in SRA Toolkit.

```bash
cd /home/data \
    && fastq-dump --verbose --split-files ERR037801
```

and also the drug resistance genes set (MEGARes, doi:10.1093/nar/gkw1009).

```bash
cd ../MEGAres \
    && wget https://megares.meglab.org/download/megares_v1.01/megares_database_v1.01.fasta
```

### Filtering

First perform the filtering of data using [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) Java-based software.

Parameters:

* `PE`: for paired-end sequencing files;

* `-phred33`: indicates the quality scoring scheme;

* `ERR037801_*.fastq`: the input files;

* `ERR037801_\*\_FILTERED.fastq` and `ERR037801_\*\_UNPAIRED.fastq`: paired and unpaired of forward and reverse reads respectively;

* `ILLUMINACLIP`: specify the adapters file. In this case the the adapters set available by Trimmomatic was used, but a set of custom adapters would be provided instead;

* `LEADING`: Remove low quality bases from the beginning. Specifies the minimum quality required to keep a base;

* `TRAILING`: Remove low quality bases from the end. Specifies the minimum quality required to keep a base;

* `SLIDINGWINDOW`: Perform a sliding window trimming, cutting once the average quality within the window falls
below a threshold. By considering multiple bases, a single poor quality base will not cause the
removal of high quality data later in the read. 

* `MINLEN`: removes reads that fall below the specified minimal length. Specifies the minimum length of reads to be kept;

* `-validatePairs`: validate if sequences keep their pairs after filtering process (maybe removed at latest versions).


```bash
java -jar \
    /home/bioinfo/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE -phred33 \
    /home/data/ERR037801_1.fastq \
    /home/data/ERR037801_2.fastq \
    /home/data/ERR037801_1_FILTERED.fastq \
    /home/data/ERR037801_1_UNPAIRED.fastq \
    /home/data/ERR037801_2_FILTERED.fastq \
    /home/data/ERR037801_2_UNPAIRED.fastq \
    ILLUMINACLIP:/home/bioinfo/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:50:20 \
    MINLEN:36 \
    -validatePairs
```

OBS: Then [read the documentation](http://www.usadellab.org/cms/?page=trimmomatic): see the introduction from Trimmomatic manual.

Trimmomatic is a fast, multi-threaded command line tool that can be used to trim and crop
Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
depending on the library preparation and downstream application.

There are two major modes of the program: Paired end mode and Single end mode. The
paired end mode will maintain correspondence of read pairs and also use the additional
information contained in paired reads to better find adapter or PCR primer fragments
introduced by the library preparation process.

Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores,
depending on the Illumina pipeline used). Files compressed using either `.gzip` or `.bzip2` are
supported, and are identified by use of `.gz` or`.bz2` file extensions. 

### Assembly

Sequences assembly can take a long time. This process is performed using `SPAdes`, a Python software. `SPAdes` is a de Bruijn graph-based assembler. The Assemble Reads with SPAdes App allows the user to assemble a genome from reads using the `SPAdes 3.13.0` assembler, which is designed for *small genomes* and *single cell sequencing*.

```bash
python3 \
    /home/bioinfo/SPAdes-3.13.0-Linux/bin/spades.py \
    --pe1-1 /home/data/ERR037801_1_FILTERED.fastq \
    --pe1-2 /home/data/ERR037801_2_FILTERED.fastq \
    -t 4 --careful --cov-cutoff auto \
    -o /home/data/ERR037801_assembly
```

### Get assembly stats

```bash
/home/bioinfo/quast-5.0.2/quast.py /home/data/ERR037801_assembly/contigs.fasta
```

### Resistance genes prediction

About Glimmer ([original description](http://ccb.jhu.edu/software/glimmer/index.shtml)): Glimmer is a system for finding genes in microbial DNA, especially the genomes of bacteria, archaea, and viruses. Glimmer (Gene Locator and Interpolated Markov ModelER) uses interpolated Markov models (IMMs) to identify the coding regions and distinguish them from noncoding DNA. The IMM approach is described in our original [Nucleic Acids Research paper on Glimmer 1.0](http://ccb.jhu.edu/papers/glimmer-nar.pdf) and in our [subsequent paper on Glimmer 2.0](http://ccb.jhu.edu/papers/glimmer2.pdf). The IMM is a combination of Markov models from 1st through 8th-order, where the order used is determined by the amount of data available to train the model. In addition, the positions used as context for the model need not immediately precede the predicted postion but are determined by a decision procedure based on the predictive power of each position in the training data set (which we term an Interpolated Context Model or ICM). The models for coding sequence are 3-periodic nonhomogenous Markov models. Improvements made in version 3 of Glimmer are described in the [third Glimmer paper](http://ccb.jhu.edu/papers/glimmer3.pdf).

First constructs the interpolated context model (ICM).

```bash
/home/bioinfo/glimmer3.02/bin/build-icm \
    MEGARes < megares_database_v1.01.fasta
```

Next, predict genes based on MEGAres binary file.

```bash
/home/bioinfo/glimmer3.02/bin/glimmer3 \
    --linear \
    ERR037801_assembly/contigs.fasta \
    MEGARes \
    drug_resistence_genes
```

And, store results in a fasta file.

```bash
/home/bioinfo/glimmer3.02/bin/extract \
    ERR037801_assembly/contigs.fasta \
    drug_resistence_genes.predict > \
    putative_drug_resistence_genes.fasta
```

### Gene annotation

Annotate contigs based in the MEGAres database. Then, first build a reference database using:

```bash
makeblastdb \
    -dbtype nucl \
    -in /home/MEGAres/megares_database_v1.01.fasta \
    -out /home/blast/MEGARes
```

and performing a blast comparison:

```bash
blastn \
    -query putative_drug_resistence_genes.fasta \
    -evalue 10e-5 -db MEGARes \
    -out putative_drug_resistence_genes_annotated.csv
```

You can optionally export the blastn results in a csv file including the `-outfmt=10` parameter (the number after equal signal can be change according the desired format, see documentation or simple use `blastn -help`).

### Extract results from container (optional)

Results from docker container are stored inside it. You can copy them to a local repository running:

`docker cp <CONTAINER_ID>:<CONTAINER_PATH> <LOCAL_DIRECTORY_PATH>`

* `docker cp`: the docker software call and the linux copy command (`cp`) to execute the action;

* `CONTAINER_ID`: the container ID recovered using `docker ps`;

* `CONTAINER_PATH`: the directory of docker image containing the files to recovery;

* `LOCAL_DIRECTORY_PATH`: the local directory to store the results files.

Example:

```bash
docker cp 90ffb38c27c9:home/ ./path/
```

A shell script containing all code explained above combined in a single pipeline is available at the `02_de_novo_assembly.sh` file located in the `~/sequence-assembly/scripts/` directory. It can be used to automated run all steps.

NOTE: some divergences can occur between this file and the shell script (`02_de_novo_assembly.sh`).
