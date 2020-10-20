#!/bin/bash
#
# NAME: De-novo assembly pipeline
# PLAT: Linux
# DESCRIPTION: Contain pipeline to perform a De-novo sequences assembly and
# Resistence genes annotation give a reference genome.
#
# ERR037801 used for tests.

# Set the timezone
export TZ=/usr/share/zoneinfo/America/Sao_Paulo

# Source base log colors
source "./colors.sh"

: "# ************************************************ #"
: "Set the base work directory."

BASE_DIR='/home'

: "Verify it the input was passed as parameter."

if [[ ! $1 ]]; then
	echo "Need at last one parameter."
	exit 1
else
	TARGET=$1
fi

echo -e "${SUCCESS} Target: $TARGET"

: "# ************************************************ #"
: "Create the control file. It should be used in flux control of pipeline."

CONTROL_FILE="${BASE_DIR}/.control-file"
source $CONTROL_FILE

[[ -z $DENOVO_DOWNLOAD_TARGET ]] && DENOVO_DOWNLOAD_TARGET=0
[[ -z $DENOVO_DOWNLOAD_REFERENCE ]] && DENOVO_DOWNLOAD_REFERENCE=0
[[ -z $DENOVO_FILTERING ]] && DENOVO_FILTERING=0
[[ -z $DENOVO_PREDICT ]] && DENOVO_PREDICT=0
[[ -z $DENOVO_ANNOTATE ]] && DENOVO_ANNOTATE=0

echo -e "\nConcluded steps (0=pending, 1=concluded):"
echo -e "\t Download targed: $DENOVO_DOWNLOAD_TARGET"
echo -e "\t Download reference: $DENOVO_DOWNLOAD_REFERENCE"
echo -e "\t Filtering: $DENOVO_FILTERING"
echo -e "\t Predict: $DENOVO_PREDICT"
echo -e "\t Annotate: $DENOVO_ANNOTATE"

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

: "# ************************************************ #"
: "Check work directories structure."

echo -e "${INFO} Check work directories structure"

MESSAGE="directory exists. Please specify a new base directory defining 'BASE_DIR' variable."

[ -d $BASE_DIR ] &&
	echo -e "${WARN} Base directory already exists." ||
	mkdir $BASE_DIR

cd $BASE_DIR

MEGARES="MEGAres"
[ -d ${MEGARES} ] &&
	echo -e "${WARN} MEGAres ${MESSAGE}" ||
	mkdir ${MEGARES}

BLAST="blast"
[ -d ${BLAST} ] &&
	echo -e "${WARN} blast ${MESSAGE}" ||
	mkdir ${BLAST}

DATA="data"
[ -d ${DATA} ] &&
	echo -e "${WARN} data ${MESSAGE}" ||
	mkdir ${DATA}

GLIMMER="glimmer"
[ -d ${GLIMMER} ] &&
	echo -e "${WARN} glimmer ${MESSAGE}" ||
	mkdir ${GLIMMER}

ls -l

: "# ************************************************ #"
: "Downloading record from Genbank."

echo -e "${INFO} Download target record"

if [[ $DENOVO_DOWNLOAD_TARGET = 0 ]]; then
	cd ${DATA} &&
		fastq-dump --split-files $TARGET &&
		cd ..

	set_control_variables 0 DENOVO_DOWNLOAD_TARGET
else
	set_control_variables 1 DENOVO_DOWNLOAD_TARGET
fi

: "# ************************************************ #"
: "Download drug resistance genes set of MEGARes (doi:10.1093/nar/gkw1009)"

echo -e "${INFO} Download reference record"

if [[ $DENOVO_DOWNLOAD_REFERENCE = 0 ]]; then
	cd ${MEGARES} &&
		wget https://megares.meglab.org/download/megares_v1.01/megares_database_v1.01.fasta &&
		cd ..

	set_control_variables 0 DENOVO_DOWNLOAD_REFERENCE
else
	set_control_variables 1 DENOVO_DOWNLOAD_REFERENCE
fi

: "# ************************************************ #"
: "Filtering and assembly."

echo -e "${INFO} Filtering, assembly, and get stats"

if [[ $DENOVO_FILTERING = 0 ]]; then
	cd ${DATA} &&
		java -jar \
			../bio-softwares/Trimmomatic-0.39/trimmomatic-0.39.jar \
			PE -phred33 \
			${TARGET}_1.fastq \
			${TARGET}_2.fastq \
			${TARGET}_1_FILTERED.fastq \
			${TARGET}_1_UNPAIRED.fastq \
			${TARGET}_2_FILTERED.fastq \
			${TARGET}_2_UNPAIRED.fastq \
			ILLUMINACLIP:../bio-softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:36 \
			-validatePairs &&
		python3 \
			../bio-softwares/SPAdes-3.13.0-Linux/bin/spades.py \
			--pe1-1 ${TARGET}_1_FILTERED.fastq \
			--pe1-2 ${TARGET}_2_FILTERED.fastq \
			-t 4 --careful --cov-cutoff auto \
			-o ${TARGET}_assembly &&
		./bio-softwares/quast-5.0.2/quast.py \
			${TARGET}_assembly/contigs.fasta &&
		cd ..

	set_control_variables 0 DENOVO_FILTERING
else
	set_control_variables 1 DENOVO_FILTERING
fi

: "# ************************************************ #"
: "Resistance genes prediction"

echo -e "${INFO} Resistance genes prediction"

if [ $DENOVO_PREDICT = 0 ]; then
	cd ${MEGARES} &&
		./bio-softwares/glimmer3.02/bin/build-icm \
			MEGARes <megares_database_v1.01.fasta &&
		./bio-softwares/glimmer3.02/bin/glimmer3 \
			--linear \
			../${DATA}/${TARGET}_assembly/contigs.fasta \
			MEGARes \
			drug_resistence_genes &&
		./bio-softwares/glimmer3.02/bin/extract \
			../${DATA}/${TARGET}_assembly/contigs.fasta \
			drug_resistence_genes.predict > \
			putative_drug_resistence_genes.fasta &&
		cd ..

	set_control_variables 0 DENOVO_PREDICT
else
	set_control_variables 1 DENOVO_PREDICT
fi

: "# ************************************************ #"
: "Gene annotation."

echo -e "${INFO} Gene annotation"

if [ $DENOVO_ANNOTATE = 0 ]; then
	cd ${BLAST} &&
		makeblastdb \
			-dbtype nucl \
			-in ../${MEGARES}/megares_database_v1.01.fasta \
			-out MEGARes &&
		blastn \
			-query putative_drug_resistence_genes.fasta \
			-evalue 10e-5 -db MEGARes \
			-out putative_drug_resistence_genes_annotated.csv &&
		cd ..

	set_control_variables 0 DENOVO_ANNOTATE
else
	set_control_variables 1 DENOVO_ANNOTATE
fi
