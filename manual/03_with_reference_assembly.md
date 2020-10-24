# Genome assembly with reference genome (NOT FINISHED)

Aim: perform a genome assembly with reference sequence.

This tutorial explains how to:

1. Get the genome of *Enterobacter cloacae* subsp. *cloacae* (NCTC 9394 isolate) from Genbank;

2. Get the MEGAres database used for the prediction of resistance genes;

3. Cleaning the SRA file (filtering step) using Trimmomatic;

4. Get a reference genome;

5. Indexing reference genome and map reads to reference;

6. And finally, generate SAM and BAM files with mapped results.

NOTE: Use the same environment created at the last tutorial.

## Downloading Sequences

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
    && wget https://megares.meglab.org/download/megares_v1.01/megares_database_v1.01.fasta \
```

## Filtering

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
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -validatePairs
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

## Assembly

### Download reference genome

First create a dedicated directory for the reference genome.

```bash
mkdir reference_genome
```

And download the reference genomes inside itusing `ncbi-genome-download` software.

Parameters:

* `bacteria`: The first parameter specify the NCBI taxonomic groups to download;
    * Options are: 'all', 'archaea', 'bacteria', 'fungi', 'invertebrate', 'metagenomes', 'plant', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral'.
    * Usage for more then one group: "bacteria,viral".

* `-g` (also --genera): specify the target genera;

* `-T` (also --species-taxids): specify the standard taxids from Genbank taxonomy;

* `-l` (also --assembly-levels): specify the assembly levels of genomes to download;
    * Options are: 'all', 'complete', 'chromosome', 'scaffold', 'contig'.
    * Usage for more then one group: "complete,chromosome".

* `-F` (also --formats): Which formats to download.
    * Options are: 'genbank', 'fasta', 'rm', 'features', 'gff', 'protein-fasta', 'genpept', 'wgs', 'cds-fasta', 'rna-fna', 'rna-fasta', 'assembly-report', 'assembly-stats', 'all'.
    * Usage for more then one group: "fasta,assembly-report".

* `-N` (also --no-cache): Don't cache the assembly summary file;

* `-v` (also --verbose): Increase output verbosity.

The full set of options can be viewed using `ncbi-genome-download --help`.

```bash
ncbi-genome-download \
    bacteria \
    -g Enterobacter \
    -T 550 \
    -l complete \
    -F fasta \
    -N -v
```

```bash
gunzip \
    refseq/bacteria/GCF_000025565.1/GCF_000025565.1_ASM2556v1_genomic.fna.gz \
    && cp refseq/bacteria/GCF_000025565.1/GCF_000025565.1_ASM2556v1_genomic.fna .
```

### Generate genome index for Burrows-Wheeler Alignment (BWA)

**Official description**: HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome, and exome sequencing data) against the general human population (as well as against a single reference genome). Based on [GCSA](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6698337&tag=1) (an extension of [BWT](http://en.wikipedia.org/wiki/Burrows-Wheeler_transform) for a graph), we designed and implemented a graph FM index (GFM), an original approach and its first implementation to the best of our knowledge. In addition to using one global GFM index that represents general population, HISAT2 uses a large set of small GFM indexes that collectively cover the whole genome (each index representing a genomic region of 56 Kbp, with 55,000 indexes needed to cover human population). These small indexes (called local indexes) combined with several alignment strategies enable effective alignment of sequencing reads. This new indexing scheme is called Hierarchical Graph FM index (HGFM). We have developed HISAT 2 based on the [HISAT](http://ccb.jhu.edu/software/hisat) and [Bowtie2](http://bowtie-bio.sf.net/bowtie2) implementations. HISAT2 outputs alignments in [SAM](http://samtools.sourceforge.net/SAM1.pdf) format, enabling interoperation with a large number of other tools (e.g. [SAMtools](http://samtools.sourceforge.net/), [GATK](http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit)) that use SAM. HISAT2 is distributed under the [GPLv3 license](http://www.gnu.org/licenses/gpl-3.0.html), and it runs on the command line under Linux, Mac OS X and Windows.

```bash
hisat2-build \
    GCF_000025565.1_ASM2556v1_genomic.fna \
    GCF_000025565.1_ASM2556v1_genomic.idx
```

### Map reads to the reference genome

Read about the output format ([SAM](https://en.wikipedia.org/wiki/SAM_(file_format)))

```bash
hisat2 \
    -1 ERR037801_1_FILTERED.fastq \
    -2 ERR037801_2_FILTERED.fastq \
    --summary-file summary.txt -p 4 \
    --no-discordant --no-mixed \
    -x GCF_000025565.1_ASM2556v1_genomic.idx \
    -S hisat2_result.sam
```

### Convert sam to bam and sort

```bash
samtools view \
    -bS hisat2_result.sam | samtools \
    sort -o hisat2_result.bam
```

### Create index for the sorted bam

```bash
samtools index \
    hisat2_result.bam \
    hisat2_result.bai
```

### Count mapped reads

```bash
samtools view \
    -bS hisat2_result.sam \
    | cut -f1 \
    | sort \
    | uniq \
    | wc -l
```

### SVN

```bash
samtools mpileup \
    -f GCF_000025565.1_ASM2556v1_genomic.fna \
    -v -u \
    hisat2_result.bam > hisat2_result.mpileup
```

```bash
sed -n '500,515p' hisat2_result.mpileup
```

Example:

```bash
docker cp 90ffb38c27c9:home/ ./path/
```

A shell script containing all code explained above combined in a single pipeline is available at the `03_with_reference_assembly.sh` file located in the `~/sequence-assembly/scripts/` directory. It can be used to automated run all steps.

NOTE: some divergences can occur between this file and the shell script (`03_with_reference_assembly.sh`).
