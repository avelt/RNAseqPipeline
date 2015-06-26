#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# This script creates a graph showing the coverage of reads along the genes with RSeQC
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh GeneBodyCoverage.sh -c config-file -o outdir [-h]

		Description:
		    Creation of a graph showing the coverage of reads along the genes with RSeQC

		Options:
			-c, --config-file Config file of the run.
		        This option is required. 
		        Generally named : config.sh
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

echo "****************************************************************"
echo "`date`: generation of graphic of gene body coverage : start"
echo "****************************************************************"

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

cd ${BIN}/${RSEQC_VERSION}
export PYTHONPATH=${PYTHONPATH}/${RSEQC_PYTHONPATH}/

bedfile=$( create_bed_from_gtf ${GTF_FILE} ${SCRIPTSDIR} )

# list containing all the bam files
ls ${OUTDIR}/Alignment/*/*_alignment.bam > ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/list_of_bam_files.txt

# if we use furious, we use softwareSL and python 2.7 must be recovered (else python 2.6 is used bu default)
echo $BIN | if grep -q "softwareSL\|biopuces"; then
	
	echo "${BIN}/Python/bin/python geneBody_coverage.py -r ${bedfile} -i ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/list_of_bam_files.txt -o ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage"
	${BIN}/Python/bin/python geneBody_coverage.py -r ${bedfile} -i ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/list_of_bam_files.txt -o ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage

else

	echo "python geneBody_coverage.py -r ${bedfile} -i ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/list_of_bam_files.txt -o ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage"
	python geneBody_coverage.py -r ${bedfile} -i ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/list_of_bam_files.txt -o ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage

fi

cd ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/
mv ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage*.r ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage_plot.r
wait
sed -i "s/legend(0,1,/legend(0,1,cex=0.6,/g" ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage_plot.r
sed -i "s/scale=c(\"none\"),/scale=c(\"none\"),cex=0.6,/g" ${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage_plot.r

echo "source(\"${OUTDIR}/Quality/RSeQC/geneBodyCoverage/geneBodyCoverage_plot.r\")" | $RBIN --no-save --no-restore --quiet

echo "****************************************"
echo "`date`: generation of graphic of gene body coverage : done"
echo "****************************************"
