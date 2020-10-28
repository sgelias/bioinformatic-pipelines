#!/bin/bash
#
# NAME: Salmonella serovar identification pipeline
# PLAT: Linux
# DESCRIPTION: Contain pipeline to perform Salmonella in silico
# serovar identification.
#
# Based in https://doi.org/10.3389/fmicb.2019.00835.

# Set the timezone
export TZ=/usr/share/zoneinfo/America/Sao_Paulo

# Source base log colors
source "./colors.sh"

: "# ************************************************ #"
: "Verify it the input was passed as parameter."

while test $# -gt 0; do
    case "$1" in
    -h | --help)
        echo "$package - attempt to capture frames"
        echo " "
        echo "$package [options] application [arguments]"
        echo " "
        echo "options:"
        echo "-h, --help                show brief help"
        echo "-a, --action=ACTION       specify an action to use"
        echo "-o, --output-dir=DIR      specify a directory to store output in"
        exit 0
        ;;
    -a)
        shift
        if test $# -gt 0; then
            export PROCESS=$1
        else
            echo "no process specified"
            exit 1
        fi
        shift
        ;;
    --action*)
        export PROCESS=$(echo $1 | sed -e 's/^[^=]*=//g')
        shift
        ;;
    -o)
        shift
        if test $# -gt 0; then
            export OUTPUT=$1
        else
            echo "no output dir specified"
            exit 1
        fi
        shift
        ;;
    --output-dir*)
        export OUTPUT=$(echo $1 | sed -e 's/^[^=]*=//g')
        shift
        ;;
    *)
        break
        ;;
    esac
done



if [[ ! $1 ]]; then
	echo "Need at last one parameter."
	exit 1
else
	TARGET=$1
fi

if [[ ! $1 ]]; then
	echo "Need at last one parameter."
	exit 1
else
	TARGET=$1
fi

echo -e "${SUCCESS} Target: $TARGET"

: "# ************************************************ #"
: "Set base variables."

BASE_DIR='/home/serovar'
SOFTWARES_DIR='/home/bio-softwares'
DATA_DIR='/home/data/salmonella-serovar'

: "# ************************************************ #"
: "Create the control file. It should be used in flux control of pipeline."

CONTROL_FILE="${BASE_DIR}/.control-file"

if [ ! -a $CONTROL_FILE ]; then
    touch $CONTROL_FILE
fi

source $CONTROL_FILE

[[ -z $SEROVAR_DOWNLOAD_TARGET ]] && SEROVAR_DOWNLOAD_TARGET=0
[[ -z $SEROVAR_ASSEMBLY ]] && SEROVAR_ASSEMBLY=0
[[ -z $SEROVAR_BUILD_SEROVAR_DATABASE ]] && SEROVAR_BUILD_SEROVAR_DATABASE=0

echo -e "\nConcluded steps (0=pending, 1=concluded):"
echo -e "\t Download targed: $SEROVAR_DOWNLOAD_TARGET"
echo -e "\t Assembly: $SEROVAR_ASSEMBLY"
echo -e "\t Build serovar datasbase: $SEROVAR_BUILD_SEROVAR_DATABASE"

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

batch_download_sra() {
    # Given a list of SRA accessions get each one from Genbank.

    CONCLUDED_DOWNLOADS="./.concluded"
    SEROVAR_FILE=$1

    if [ ! -a $CONCLUDED_DOWNLOADS ]; then
        touch $CONCLUDED_DOWNLOADS
    fi

    while read line; do
        if [[ -z $(grep $line $CONCLUDED_DOWNLOADS) ]]; then
            echo "Downloading $line ..."
            fastq-dump --split-files --fasta "default" $line
            echo $line >>$CONCLUDED_DOWNLOADS
        else
            echo "$line record also downloaded. Ignored."
        fi
    done <$SEROVAR_FILE
}

populate_reference_dataset() {
    # Populate reference fasta files into a single file.

    FILEPATH="${DATA_DIR}/serovar-specific-sequences"
    SEROVAR_FILE=$1

    if [ ! -d $FILEPATH ]; then
        echo "Invalid reference dataset path."
        exit 1
    fi

    if [ -a $SEROVAR_FILE ]; then
        >$SEROVAR_FILE
    fi

    for filename in $(find ${FILEPATH} -type f -name "*.fasta"); do
        [ -e "$filename" ] || continue
        cat $filename >>$SEROVAR_FILE
    done
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
BLAST="blast"

[ -d ${ASSEMBLY} ] &&
    echo -e "${WARN} data ${MESSAGE}" ||
    mkdir ${ASSEMBLY}

[ -d ${BLAST} ] &&
    echo -e "${WARN} blast ${MESSAGE}" ||
    mkdir ${BLAST}

: "# ************************************************ #"
: "Downloading record from Genbank."

echo -e "${INFO} Download target record"

if [[ $SEROVAR_DOWNLOAD_TARGET = 0 ]]; then
    cd ${DATA_DIR} &&
        batch_download_sra $TARGET &&
        cd ${BASE_DIR}

    set_control_variables 0 SEROVAR_DOWNLOAD_TARGET
else
    set_control_variables 1 SEROVAR_DOWNLOAD_TARGET
fi

: "# ************************************************ #"
: "Assembly."



# TERMINAR DE FAZER O ASSEMBLY COM SPADES



echo -e "${INFO} Assembly, and get stats"

if [[ $SEROVAR_ASSEMBLY = 0 ]]; then
	cd ${BASE_DIR}/${ASSEMBLY} &&
		python3 \
			${SOFTWARES_DIR}/SPAdes-3.13.0-Linux/bin/spades.py \
			--pe1-1 ${TARGET}_1.fastq \
			--pe1-2 ${TARGET}_1.fastq \
			-t 4 --careful --cov-cutoff auto \
			-o ${TARGET}_assembly &&
		${SOFTWARES_DIR}/quast-5.0.2/quast.py \
			${TARGET}_assembly/contigs.fasta &&
		cd ${BASE_DIR}

	set_control_variables 0 SEROVAR_ASSEMBLY
else
	set_control_variables 1 SEROVAR_ASSEMBLY
fi

: "# ************************************************ #"
: "Build serovar specific blast database."

echo -e "${INFO} Build serovar specific blast database"

if [[ $SEROVAR_BUILD_SEROVAR_DATABASE = 0 ]]; then
    SEROVAR_FILE="./serovar.fasta"

    cd ${BASE_DIR}/${BLAST} &&
        populate_reference_dataset ${SEROVAR_FILE} &&
        makeblastdb \
            -dbtype nucl \
            -in ${SEROVAR_FILE} \
            -out serovar &&
        cd ${BASE_DIR}

    set_control_variables 0 SEROVAR_BUILD_SEROVAR_DATABASE
else
    set_control_variables 1 SEROVAR_BUILD_SEROVAR_DATABASE
fi
