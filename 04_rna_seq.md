# Explore RNAseq results

Aim: explore RNAseq results.

This tutorial explains how to:

1. Pull a docker environment containing all tools necessary to perform the genome assembly and annotation;

2. Prepare data;

3. Sequences assembly;

4. Validate assembled genes using a reference organism (*Saccharomyces pombe*);

5. Get statistics of the gene expression.

## Pull docker image using pure docker or docker-compose (optional)

The docker image should be pulled from Docker Hub using the command below. Alternatively, a non-docker Linux machine can be prepared following the steps described in the section Environment preparation.

```bash
docker pull waldeyr/rnaseq:1.0
```

Then, turn on the docker container using:

```bash
docker run --rm -it waldeyr/rnaseq:1.0
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
        container_name: "rna_seq"
        image: waldeyr/rnaseq:1.0
        volumes:
            - analysis:/home/
```

Then, start the container using:

```bash
docker-compose run --rm machine
```

Using the docker-compose you have several advantages, including persist all parameters or/and create once all containers that composes the *orchestra*.

## Build analysis data

Concatenate data

```bash
cat *.left.fq > reads.ALL.left.fq && cat *.right.fq > reads.ALL.right.fq
```

## Running analyses

### Assembly and assembly stats

```bash
Trinity \
    --CPU 4 --max_memory 6G \
    --seqType fq \
    --SS_lib_type RF \
    --left reads.ALL.left.fq \
    --right reads.ALL.right.fq
```

```bash
$TRINITY_HOME/util/TrinityStats.pl \
    trinity_out_dir/Trinity.fasta
```

### Checking transcripts matches among this study and the reference organism with Blast

```bash
makeblastdb \
    -in GCF_000002945.1_ASM294v2_rna.fna \
    -dbtype nucl
```

```bash
blastn \
    -query trinity_out_dir/Trinity.fasta \
    -db GCF_000002945.1_ASM294v2_rna.fna \
    -out Trinity_vs_S_pombe_rna.blastn \
    -evalue 1e-20 \
    -dust no \
    -task megablast \
    -num_threads 4 \
    -max_target_seqs 1 \
    -outfmt 6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl \
    Trinity_vs_S_pombe_rna.blastn \
    trinity_out_dir/Trinity.fasta \
    GCF_000002945.1_ASM294v2_rna.fna
```

### Abundance estimation

```bash
$TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcripts trinity_out_dir/Trinity.fasta \
    --seqType fq \
    --left Sp.ds.1M.left.fq \
    --right Sp.ds.1M.right.fq \
    --est_method RSEM \
    --aln_method bowtie \
    --trinity_mode \
    --prep_reference \
    --output_dir rsem_DS_outdir
```

```bash
$TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcripts trinity_out_dir/Trinity.fasta \
    --seqType fq \
    --left Sp.hs.1M.left.fq \
    --right Sp.hs.1M.right.fq \
    --est_method RSEM \
    --aln_method bowtie \
    --trinity_mode \
    --prep_reference \
    --output_dir rsem_HS_outdir
```

```bash
$TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcripts trinity_out_dir/Trinity.fasta \
    --seqType fq --left Sp.log.1M.left.fq \
    --right Sp.log.1M.right.fq \
    --est_method RSEM \
    --aln_method bowtie \
    --trinity_mode \
    --prep_reference \
    --output_dir rsem_LOG_outdir
```

```bash
$TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcripts trinity_out_dir/Trinity.fasta \
    --seqType fq \
    --left Sp.plat.1M.left.fq \
    --right Sp.plat.1M.right.fq \
    --est_method RSEM \
    --aln_method bowtie \
    --trinity_mode \
    --prep_reference \
    --output_dir rsem_PLAT_outdir
```

#### Matrix

```bash
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map \
    --est_method RSEM \
    --name_sample_by_basedir \
    rsem_DS_outdir/RSEM.isoforms.results \
    rsem_HS_outdir/RSEM.isoforms.results \
    rsem_LOG_outdir/RSEM.isoforms.results \
    rsem_PLAT_outdir/RSEM.isoforms.results
```

#### Plot Heatmap

```bash
Rscript plotHeatMap.R
```
