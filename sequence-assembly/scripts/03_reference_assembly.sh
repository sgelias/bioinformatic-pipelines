#!/bin/bash
#
# NAME: Reference assembly pipeline
# PLAT: Linux
# DESCRIPTION: Contain pipeline to perform a sequences assembly and
# using a reference genome and annotate resistance genes give a
# reference genome.
#
# ERR037801 used for tests.

# Set the timezone
export TZ=/usr/share/zoneinfo/America/Sao_Paulo

# Source base log colors
source "./colors.sh"

: "# ************************************************ #"
: "Verify it the input was passed as parameter."

if [[ ! $1 ]]; then
    echo "Need at last one parameter."
    exit 1
else
    TARGET=$1
fi

echo -e "${SUCCESS} Target: $TARGET"

: "# ************************************************ #"
: "Set base variables."

BASE_DIR='/home/ref-map'
SOFTWARES_DIR='/home/bio-softwares'
DATA_DIR='/home/data'

: "# ************************************************ #"
: "Create the control file. It should be used in flux control of pipeline."

CONTROL_FILE="${BASE_DIR}/.control-file"
source $CONTROL_FILE

[[ -z $REFMAP_DOWNLOAD_TARGET ]] && REFMAP_DOWNLOAD_TARGET=0
[[ -z $REFMAP_DOWNLOAD_REFERENCE ]] && REFMAP_DOWNLOAD_REFERENCE=0
[[ -z $REFMAP_FILTERING ]] && REFMAP_FILTERING=0
[[ -z $REFMAP_MAPTOREFERENCE ]] && REFMAP_MAPTOREFERENCE=0
#[[ -z $REFMAP_ANNOTATE ]] && REFMAP_ANNOTATE=0

echo -e "\nConcluded steps (0=pending, 1=concluded):"
echo -e "\t Download targed: $REFMAP_DOWNLOAD_TARGET"
echo -e "\t Download reference: $REFMAP_DOWNLOAD_REFERENCE"
echo -e "\t Filtering: $REFMAP_FILTERING"
echo -e "\t Predict: $REFMAP_MAPTOREFERENCE"
#echo -e "\t Annotate: $REFMAP_ANNOTATE"

: "# ************************************************ #"
: "Create auxiliary functions"

set_control_variables() {
    # Called when control flux variables would be changed.
    # First argument should be 0 if the current step was resolved
    # in the current run 1 if already resolved in previous run.

    VAR_DATE="${2}_DATE"

    if [ "$1" = 0 ]; then
        echo "${2}=1" >>$CONTROL_FILE
        echo "${VAR_DATE}=$(date +%F:%H-%M-%S:%Z)" >>$CONTROL_FILE
        echo -e "${SUCCESS} Target download already concluded in: $(date +%F:%H-%M-%S:%Z)"
    elif [ "$1" = 1 ]; then
        echo -e "${SUCCESS} Target download already concluded in: ${!VAR_DATE}"
    fi
}

: "Check work directories structure."

echo -e "${INFO} Check work directories structure"

MESSAGE="directory exists. Please specify a new base directory defining 'BASE_DIR' variable."

[ -d $BASE_DIR ] &&
    echo -e "${WARN} Base directory already exists." ||
    mkdir $BASE_DIR

cd $BASE_DIR

ASSEMBLY="assembly"
REFERENCE="reference"
MAPPED="mapped"

[ -d ${ASSEMBLY} ] &&
    echo -e "${WARN} data ${MESSAGE}" ||
    mkdir ${ASSEMBLY}

[ -d ${REFERENCE} ] &&
    echo -e "${WARN} reference ${MESSAGE}" ||
    mkdir ${REFERENCE}

[ -d ${MAPPED} ] &&
    echo -e "${WARN} mapped ${MESSAGE}" ||
    mkdir ${MAPPED}

: "# ************************************************ #"
: "Downloading record from Genbank."

echo -e "${INFO} Download target record"

if [[ $REFMAP_DOWNLOAD_TARGET = 0 ]]; then
    cd ${DATA_DIR} &&
        fastq-dump --split-files $TARGET &&
        cd ${BASE_DIR}

    set_control_variables 0 REFMAP_DOWNLOAD_TARGET
else
    set_control_variables 1 REFMAP_DOWNLOAD_TARGET
fi

: "# ************************************************ #"
: "Download available complete genomes of Enterobacter from Genbank"

echo -e "${INFO} Download reference genome"

if [[ $REFMAP_DOWNLOAD_REFERENCE = 0 ]]; then
    cd ${BASE_DIR}/${REFERENCE} &&
        ncbi-genome-download \
            bacteria \
            -g Enterobacter \
            -T 550 \
            -l complete \
            -F fasta \
            -N -v &&
        gunzip \
            refseq/bacteria/GCF_000025565.1/GCF_000025565.1_ASM2556v1_genomic.fna.gz &&
        cp refseq/bacteria/GCF_000025565.1/GCF_000025565.1_ASM2556v1_genomic.fna . &&
        cd ${BASE_DIR}

    set_control_variables 0 REFMAP_DOWNLOAD_REFERENCE
else
    set_control_variables 1 REFMAP_DOWNLOAD_REFERENCE
fi

: "# ************************************************ #"
: "Filtering and assembly."

echo -e "${INFO} Filtering, assembly, and get stats"

if [[ $REFMAP_FILTERING = 0 ]]; then
    cd ${BASE_DIR}/${ASSEMBLY} &&
        java -jar \
            ${SOFTWARES_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar \
            PE -phred33 \
            -threads 6 \
            ${DATA_DIR}/${TARGET}_1.fastq \
            ${DATA_DIR}/${TARGET}_2.fastq \
            ${TARGET}_1_FILTERED.fastq \
            ${TARGET}_1_UNPAIRED.fastq \
            ${TARGET}_2_FILTERED.fastq \
            ${TARGET}_2_UNPAIRED.fastq \
            ILLUMINACLIP:${SOFTWARES_DIR}/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            -validatePairs &&
        python3 \
            ${SOFTWARES_DIR}/SPAdes-3.13.0-Linux/bin/spades.py \
            --pe1-1 ${TARGET}_1_FILTERED.fastq \
            --pe1-2 ${TARGET}_2_FILTERED.fastq \
            -t 4 --careful --cov-cutoff auto \
            -o ${TARGET}_assembly &&
        ${SOFTWARES_DIR}/quast-5.0.2/quast.py \
            ${TARGET}_assembly/contigs.fasta &&
        cd ${BASE_DIR}

    set_control_variables 0 REFMAP_FILTERING
else
    set_control_variables 1 REFMAP_FILTERING
fi

: "# ************************************************ #"
: "Map to reference."

echo -e "${INFO} Map to reference"

if [[ $REFMAP_MAPTOREFERENCE = 0 ]]; then
    cd ${BASE_DIR}/${MAPPED} &&
        hisat2-build \
            GCF_000025565.1_ASM2556v1_genomic.fna \
            GCF_000025565.1_ASM2556v1_genomic.idx &&
        hisat2 \
            -1 ERR037801_1_FILTERED.fastq \
            -2 ERR037801_2_FILTERED.fastq \
            --summary-file summary.txt -p 4 \
            --no-discordant --no-mixed \
            -x GCF_000025565.1_ASM2556v1_genomic.idx \
            -S hisat2_result.sam &&
        samtools view \
            -bS hisat2_result.sam | samtools \
            sort -o hisat2_result.bam &&
        samtools index \
            hisat2_result.bam \
            hisat2_result.bai &&
        samtools view \
            -bS hisat2_result.sam | cut -f1 | sort | uniq | wc -l &&
        cd ${BASE_DIR}

    set_control_variables 0 REFMAP_MAPTOREFERENCE
else
    set_control_variables 1 REFMAP_MAPTOREFERENCE
fi
