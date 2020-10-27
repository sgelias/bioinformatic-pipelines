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
: "Verify it the input was passed as parameter."

if [ $# -eq 0 ]; then
	echo "No parameter provided."
	exit 1
fi

while test $# -gt 0; do
	case "$1" in
	-h | --help)
		echo "$package - attempt to capture frames"
		echo " "
		echo "$package [options] application [arguments]"
		echo " "
		echo "options:"
		echo "-h, --help                show brief help"
		echo "-t, --target=SRA       specify the target SRA to use"
		#echo "-o, --output-dir=DIR      specify a directory to store output in"
		exit 0
		;;
	-t)
		shift
		if test $# -gt 0; then
			export TARGET=$1
		else
			echo "no target specified"
			exit 1
		fi
		shift
		;;
	--target*)
		export TARGET=$(echo $1 | sed -e 's/^[^=]*=//g')
		shift
		;;
	*)
		break
		;;
	esac
done

if [ -z $TARGET ]; then
	echo -e "${WARN} Please specify the target parameter using -t or --target flag.
	Example: -t ERR037801\n"
	exit 1
else
	echo -e "${SUCCESS} Target: $TARGET"
fi

exit 1

: "# ************************************************ #"
: "Set base variables."

BASE_DIR='/home/de-novo'
SOFTWARES_DIR='/home/bio-softwares'
DATA_DIR='/home/data/enterobacter'

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
		echo -e "${SUCCESS} Step already concluded in: $(date +%F:%H-%M-%S:%Z)"
	elif [ "$1" = 1 ]; then
		echo -e "${SUCCESS} Step already concluded in: ${!VAR_DATE}"
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

ASSEMBLY="assembly"
REFERENCE="reference"
GLIMMER="glimmer"
BLAST="blast"

[ -d ${ASSEMBLY} ] &&
	echo -e "${WARN} data ${MESSAGE}" ||
	mkdir ${ASSEMBLY}

[ -d ${REFERENCE} ] &&
	echo -e "${WARN} reference ${MESSAGE}" ||
	mkdir ${REFERENCE}

[ -d ${GLIMMER} ] &&
	echo -e "${WARN} glimmer ${MESSAGE}" ||
	mkdir ${GLIMMER}

[ -d ${BLAST} ] &&
	echo -e "${WARN} blast ${MESSAGE}" ||
	mkdir ${BLAST}

: "# ************************************************ #"
: "Downloading record from Genbank."

echo -e "${INFO} Download target record"

if [[ $DENOVO_DOWNLOAD_TARGET = 0 ]]; then
	cd ${DATA_DIR} &&
		fastq-dump --split-files $TARGET &&
		cd ${BASE_DIR}

	set_control_variables 0 DENOVO_DOWNLOAD_TARGET
else
	set_control_variables 1 DENOVO_DOWNLOAD_TARGET
fi

: "# ************************************************ #"
: "Download drug resistance genes set of MEGARes (doi:10.1093/nar/gkw1009)"

echo -e "${INFO} Download reference record"

if [[ $DENOVO_DOWNLOAD_REFERENCE = 0 ]]; then
	cd ${BASE_DIR}/${REFERENCE} &&
		wget https://megares.meglab.org/download/megares_v1.01/megares_database_v1.01.fasta &&
		cd ${BASE_DIR}

	set_control_variables 0 DENOVO_DOWNLOAD_REFERENCE
else
	set_control_variables 1 DENOVO_DOWNLOAD_REFERENCE
fi

: "# ************************************************ #"
: "Filtering and assembly."

echo -e "${INFO} Filtering, assembly, and get stats"

if [[ $DENOVO_FILTERING = 0 ]]; then
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

	set_control_variables 0 DENOVO_FILTERING
else
	set_control_variables 1 DENOVO_FILTERING
fi

: "# ************************************************ #"
: "Resistance genes prediction"

echo -e "${INFO} Resistance genes prediction"

if [ $DENOVO_PREDICT = 0 ]; then
	cd ${BASE_DIR}/${REFERENCE} &&
		${SOFTWARES_DIR}/glimmer3.02/bin/build-icm \
			MEGARes <megares_database_v1.01.fasta &&
		${SOFTWARES_DIR}/glimmer3.02/bin/glimmer3 \
			--linear \
			${BASE_DIR}/${ASSEMBLY}/${TARGET}_assembly/contigs.fasta \
			MEGARes \
			drug_resistence_genes &&
		${SOFTWARES_DIR}/glimmer3.02/bin/extract \
			${BASE_DIR}/${ASSEMBLY}/${TARGET}_assembly/contigs.fasta \
			drug_resistence_genes.predict > \
			putative_drug_resistence_genes.fasta &&
		cd ${BASE_DIR}

	set_control_variables 0 DENOVO_PREDICT
else
	set_control_variables 1 DENOVO_PREDICT
fi

: "# ************************************************ #"
: "Gene annotation."

echo -e "${INFO} Gene annotation"

if [ $DENOVO_ANNOTATE = 0 ]; then
	cd ${BASE_DIR}/${BLAST} &&
		makeblastdb \
			-dbtype nucl \
			-in ${BASE_DIR}/${REFERENCE}/megares_database_v1.01.fasta \
			-out MEGARes &&
		blastn \
			-query ${BASE_DIR}/${REFERENCE}/putative_drug_resistence_genes.fasta \
			-evalue 10e-5 \
			-outfmt "10" \
			-db MEGARes \
			-out putative_drug_resistence_genes_annotated &&
		cd ${BASE_DIR}

	set_control_variables 0 DENOVO_ANNOTATE
else
	set_control_variables 1 DENOVO_ANNOTATE
fi
