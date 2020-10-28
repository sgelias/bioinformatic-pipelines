#!/bin/bash
#
# NAME: Metabarcoding processing pipeline
# PLAT: Linux
# DESCRIPTION: Contain pipeline to perform a metabarcoding processing.
#
# ERR037801 used for tests.

# Set the timezone
export TZ=/usr/share/zoneinfo/America/Sao_Paulo

# Source base log colors
source "./colors.sh"

: "# ********************************** #"
: <<'STEPS'
	------- Part 0 -------

	Step: DOWNLOAD METAGENOME DATA TO PROCESSING
	Tool: SRA-Toolkit
	Description: A custom function should be written using shell or python to download files from SRA.

	------- Part 1 -------

	Step: JOIN PAIRED READS
	Tool: PEAR
	Description: Join forward and reverse sequences into a single read.

	Step: QUALITY FILTER
	Tool: FASTX-Toolkit
	Description: 

	Step: RE-INDEX REMAINING SEQUENCES AND MERGE ALL SEQUENCES INTO A SINGLE FILE
	Tool: In house Python scripts
	Description: 

	------- Part 2 -------

	Step: DEREPLICATE
	Tool: VSEARCH or a partial script from PIPITS package
	Description: 

	Step: EXTRACT ITS AND RE-ORIENTATE
	Tool: ITSx
	Description: 

	Step: RE-INFLATE
	Tool: VSEARCH or a partial script from PIPITS package
	Description: 

	------- Part 3 -------

	Step: DEREPLICATE
	Tool: VSEARCH OR A PARTIAL SCRIPTS OF PIPITS PACKAGE
	Description: 

	Step: REMOVE SINGLETONS
	Tool: VSEARCH
	Description: 

	Step: FIND OTU'S
	Tool: VSEARCH
	Description: 

	Step: REMOVE CHIMERAS
	Tool: VSEARCH
	Description: 

	Step: MAP READS ONTO OTU'S
	Tool: VSEARCH
	Description: 

	Step: ASSIGN TAXONOMY
	Tool: RDP
	Description: 

	Step: GENERATE ABUNDANCE TABLE
	Tool: In house Python scripts
	Description: 

	------- Part 4 -------

	Step: GENERATE ECOLOGICAL ANALYSIS USING R SCRIPT
	Tool: Many packages and in house R scripts
	Description: 
STEPS

: "# ************************************************ #"
: "Validate input parameters."

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
		echo "-h, --help                   show brief help"
		echo "-t, --target-file=SRA        specify a file containing target SRA accessions to use"
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
	--target-file*)
		export TARGET=$(echo $1 | sed -e 's/^[^=]*=//g')
		shift
		;;
	*)
		break
		;;
	esac
done

if [ -z $TARGET ]; then
	echo -e "${WARN} Please specify the target parameter using -t or --target-file flag.
	Example: -t SRA.txt\n"
	exit 1
else
	echo $TARGET
	if [ -f $TARGET ]; then
		echo -e "${SUCCESS} Target: $TARGET"
	fi
fi

: "# ************************************************ #"
: "Set base variables."

BASE_DIR='/home/metabarcoding'
SOFTWARES_DIR='/home/bio-softwares'
DATA_DIR='/home/data/fungal-its'

: "# ************************************************ #"
: "Create the control file. It should be used in flux control of pipeline."

CONTROL_FILE="${BASE_DIR}/.control-file"
source $CONTROL_FILE

[[ -z $METABARCODING_DOWNLOAD_TARGET ]] && METABARCODING_DOWNLOAD_TARGET=0

echo -e "\nConcluded steps (0=pending, 1=concluded):"
echo -e "\t Download targed: $METABARCODING_DOWNLOAD_TARGET"

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

: "# ************************************************ #"
: "Check work directories structure."

echo -e "${INFO} Check work directories structure"

MESSAGE="directory exists. Please specify a new base directory defining 'BASE_DIR' variable."

[ -d $BASE_DIR ] &&
	echo -e "${WARN} Base directory already exists." ||
	mkdir $BASE_DIR

cd $BASE_DIR

ASSEMBLY="assembly"
[ -d ${ASSEMBLY} ] &&
	echo -e "${WARN} data ${MESSAGE}" ||
	mkdir ${ASSEMBLY}

: "# ************************************************ #"
: "Downloading record from Genbank."

echo -e "${INFO} Download target record"

if [[ $METABARCODING_DOWNLOAD_TARGET = 0 ]]; then
	cd ${DATA_DIR} &&
		batch_download_sra "/home/scripts/$TARGET" &&
		cd ${BASE_DIR} &&
		set_control_variables 0 METABARCODING_DOWNLOAD_TARGET
else
	set_control_variables 1 METABARCODING_DOWNLOAD_TARGET
fi
