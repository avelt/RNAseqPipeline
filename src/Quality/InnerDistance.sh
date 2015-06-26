#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# This module is used to calculate the inner distance (or insert size) between two paired RNA reads. 
# The distance is the mRNA length between two paired fragments. We first determine the genomic (DNA) size
# between two paired reads: D_size = read2_start - read1_end, then
#    if two paired reads map to the same exon: inner distance = D_size
#    if two paired reads map to different exons:inner distance = D_size - intron_size
#    if two paired reads map non-exonic region (such as intron and intergenic region): inner distance = D_size
#    The inner_distance might be a negative value if two fragments were overlapped.
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh InnerDistance.sh -c config-file -s samplename -o outdir [-h]

		Description:
		    This module is used to calculate the inner distance (or insert size) between two paired RNA reads.

		Options:
			-c, --config-file Config file of the run.
		        This option is required. 
		        Generally named : config.sh
			-s, --samplename Samplename of the current file
		        This option is required.
		 	-o, --outdir Outdir of the run
		        This option is required. 
		        Generally of type : /path/to/Num_project_Name_author
		    -h, --help
		        Print this message and exit the program.
		__EOF__
}
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Getting parameters from the input
#------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "c:s:o:h" --long "config-file:,samplename:,outdir:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-s|--samplename)
			SAMPLENAME="$2"
			shift 2 
			;;
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
[ "$SAMPLENAME" == "" ] || [ "$CONFIG_FILE" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --samplename, --config-file and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------

echo "****************************************************************"
echo "`date`: inner distance calculation of paired-end data : start"
echo "****************************************************************"

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

bedfile=$( create_bed_from_gtf ${GTF_FILE} ${SCRIPTSDIR} )

cd ${BIN}/${RSEQC_VERSION}
export PYTHONPATH=${PYTHONPATH}/${RSEQC_PYTHONPATH}/

echo "./inner_distance.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -o $OUTDIR/Quality/RSeQC/innerDistance/${SAMPLENAME} -r ${bedfile}"

# if we use furious, we use softwareSL and python 2.7 must be recovered (else python 2.6 is used bu default)
echo $BIN | if grep -q "softwareSL\|biopuces"; then

	${BIN}/Python/bin/python inner_distance.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -o $OUTDIR/Quality/RSeQC/innerDistance/${SAMPLENAME} -r ${bedfile}

else

	python inner_distance.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -o $OUTDIR/Quality/RSeQC/innerDistance/${SAMPLENAME} -r ${bedfile}

fi

cd ${OUTDIR}/Quality/RSeQC/innerDistance/
echo "source(\"${OUTDIR}/Quality/RSeQC/innerDistance/${SAMPLENAME}.inner_distance_plot.r\")" | $RBIN --no-save --no-restore --quiet

mv ${OUTDIR}/Quality/RSeQC/innerDistance/${SAMPLENAME}.inner_distance_plot.pdf ${OUTDIR}/Quality/RSeQC/innerDistance/${SAMPLENAME}_inner_distance_plot.pdf


echo "****************************************"
echo "`date`: inner distance calculation of paired-end data : done"
echo "****************************************"
