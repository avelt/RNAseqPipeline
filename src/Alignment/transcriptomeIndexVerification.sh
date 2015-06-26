#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# This script allows to check whether the transcriptome-index used by tophat exists, otherwise it is automatically created.
# Once created, it can significantly reduce the time alignment of real samples.
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh transcriptomeIndexVerification.sh -c config_file -o outdir [-h]

		Description:
		    transcriptomeIndexVerification.sh allows to check whether the transcriptome-index used by tophat2 exists,
		    otherwise it is automatically created. Once created, it can significantly reduce the time of alignment step.
		    This option is only available for Tophat2, because Tophat is not able to use the transcriptome-index.

		Options:
		 	-c, --config Config file of the run.
		        This option is required. 
		        Generally named : config.sh
		    -o, --outdir Outdir of the run (defined in the run.sh script).
		        This option is required. 
		        Generally of type : /path/to/Num_project_Name_author/Date_of_analysis
		    -h, --help
		        Print this message and exit the program.
		__EOF__
}
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Getting parameters from the input
#------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "c:o:h" --long "config:,outdir:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-c|--config)
			CONFIG_FILE="$2"
			shift 2 
			;;
		-o|--outdir)
			OUTDIR="$2"
			shift 2 
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# Checking the input parameter
#------------------------------------------------------------------------------------------------------------------
# check of config file
if ! source "$CONFIG_FILE"
then
	echo "Cannot load $CONFIG_FILE, create it and re-run transcriptomeIndexVerification.sh script." 1>&2
	exit 1 
else 
	# Getting parameters
	source "$CONFIG_FILE"
fi
# check of output directory
[ ! -d "$OUTDIR" ] && \
	echo "$OUTDIR directory doesn't exist ! Rerun with valid directory path."

source "$SCRIPTSDIR/src/Utilities/utils.sh"
#------------------------------------------------------------------------------------------------------------------

# transcriptome version is recovered via the gtf file

# Bowtie indexes path, give by the config file.
export BOWTIE_INDEXES=${BOWTIE_INDEXES}

# Version of the transcriptome, recovered from the GTF file name (so, give a suitable name to you GTF file !)
TRANSCRIPTOME_VERSION=$( basename $GTF_FILE | sed 's/.gtf//' )

if [ ! -d "${GENOME_PATH}/${SPECIES}/${GENOME_VERSION}/" ]
then
	print_err_and_exit "Cannot find ${GENOME_PATH}/${SPECIES}/${GENOME_VERSION}/, which is the path to your genome's files. Exit. Check and re-run."
fi

make_directories "${GENOME_PATH}/${SPECIES}/${GENOME_VERSION}/Tophat"

TRANSCRIPTOME_INDEX=${GENOME_PATH}/${SPECIES}/${GENOME_VERSION}/Tophat/${TRANSCRIPTOME_VERSION}/${TRANSCRIPTOME_VERSION}
TRANSCRIPTOME_PATH=`dirname $TRANSCRIPTOME_INDEX`


# path to bowtie and samtools for tophat
export PATH=$PATH:$BIN/$BOWTIE_VERSION/:$BIN/$SAMTOOLS_VERSION/

# if the transcriptome-index does not exist it is created, 
# otherwise we specifies that it exist and RNAseqPipeline continue.

if [ ! -d ${TRANSCRIPTOME_PATH} ]
then
		
	echo "****************************************" >> ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX
	echo "`date` : Creation of transcriptome-index for $TRANSCRIPTOME_VERSION version." >> ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX
	echo "****************************************" >> ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX
		
	mkdir ${TRANSCRIPTOME_PATH}
	$BIN/$TOPHAT_VERSION/tophat2 -p $NBPROC --no-coverage-search -G ${GTF_FILE} --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} >> ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX
	rm ${TRANSCRIPTOME_INDEX}.gff
	ln -s ${GTF_FILE} ${TRANSCRIPTOME_INDEX}.gff
		
	wait
		
else 
	
	echo "`date`: transcriptome-index already exist. Continue." >> ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX

fi


