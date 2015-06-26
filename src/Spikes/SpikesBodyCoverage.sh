#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# SPIKES BODY COVERAGE
# creates a graph showing the coverage of reads along spikes with RSeQC
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh SpikesBodyCoverage.sh -c config_file -o outdir [-h]

		Description:
		    This script performs all analysis on spikes, on the file ${outdir}/Alignment/${SAMPLENAME}_Spikes.sam

		Options:
		 	-c, --config-file Config file of the run.
		        This option is required. 
		        Generally named : config.sh
		 	-o, --outdir Outdir of the run (defined in the run.sh script).
		        This option is required.
		    -h, --help
		        Print this message and exit the program.
		__EOF__
}
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# Getting parameters from the input
#------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "c:o:h" --long "config-file:,outdir:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-c|--config-file)
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

#------------------------------------------------------------------------------------------------------------------
# Checking the input parameter
#------------------------------------------------------------------------------------------------------------------
[ "$CONFIG_FILE" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --config-file and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

make_directories ${OUTDIR}/Quality/Spikes/geneBodyCoverage

# list containing all the bam files
ls ${OUTDIR}/Alignment/*_without_Spikes/*_Spikes.bam > ${OUTDIR}/Quality/Spikes/geneBodyCoverage/list_of_bam_files.txt

spikesBedFile=$( create_bed_from_gtf ${GTF_SPIKES} ${SCRIPTSDIR} )

cd ${BIN}/${RSEQC_VERSION}
export PYTHONPATH=${PYTHONPATH}/${RSEQC_PYTHONPATH}

# if we use furious, we use softwareSL and python 2.7 must be recovered (else python 2.6 is used bu default)
echo $BIN | if grep -q "softwareSL\|biopuces"; then

	${BIN}/Python/bin/python geneBody_coverage.py -r ${spikesBedFile} -i ${OUTDIR}/Quality/Spikes/geneBodyCoverage/list_of_bam_files.txt -o ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples

else

	python geneBody_coverage.py -r ${spikesBedFile} -i ${OUTDIR}/Quality/Spikes/geneBodyCoverage/list_of_bam_files.txt -o ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples

fi

mv ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples*.r ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples_geneBodyCoverage_plot.r
echo "source(\"${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples_geneBodyCoverage_plot.r\")" | $RBIN --no-save --no-restore --quiet
mv ${OUTDIR}/Quality/Spikes/geneBodyCoverage/*geneBodyCoverage*pdf ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples_geneBodyCoverage_plot.pdf
