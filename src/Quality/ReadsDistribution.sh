#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# Analysis of reads distribution in the exonic , intronic and intergenic regions
# The script "read_distribution.py" generates a file with
# the number of reads in each region (intron, intergenic exonic)
# Requires the existence of the gtf reference file in bed format
# If this bed file does not exist, it is automatically created with the script gtf2bed.pl
# The script necessary to create the bed file is in the src/Utilities folder
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh ReadsDistribution.sh -c config-file -s samplename -o outdir [-h]

		Description:
		    Analysis of reads distribution in the exonic , intronic and intergenic regions

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

echo "****************************************"
echo "`date`: calculation of reads distribution : start"
echo "****************************************"

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

cd ${BIN}/${RSEQC_VERSION}
export PYTHONPATH=${PYTHONPATH}/${RSEQC_PYTHONPATH}/

readDistributionFile=${OUTDIR}/Quality/RSeQC/reads_distribution/${SAMPLENAME}_reads_distribution.txt
touch ${readDistributionFile}
resultsFileRSeQC=$OUTDIR/Quality/RSeQC/resultsTable_${SAMPLENAME}.txt

bedfile=$( create_bed_from_gtf ${GTF_FILE} ${SCRIPTSDIR} )

echo "python -u read_distribution.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -r ${bedfile} &>> ${resultsFileRSeQC}"

# if we use furious, we use softwareSL and python 2.7 must be recovered (else python 2.6 is used bu default)
echo $BIN | if grep -q "softwareSL\|biopuces"; then

	${BIN}/Python/bin/python -u read_distribution.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -r ${bedfile} &>> ${resultsFileRSeQC}

else 

	python -u read_distribution.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -r ${bedfile} &>> ${resultsFileRSeQC}

fi

wait

##########################################
# recovery of density in all the different regions (exons, introns, TSS/TES)

# recovery of density in exons (sum of density in CDS, 5'UTR and 3'UTR)
DensityExons=$( cat ${resultsFileRSeQC} | grep -iE "CDS|UTR" | cut -f4 | awk '{SUM += $4} END {print SUM}' )
# recovery of density in introns
DensityIntrons=$( cat ${resultsFileRSeQC} | grep "Introns" | cut -f3 | awk '{SUM += $4} END {print SUM}' )
# recovery of density in TSS_up_10kb and TES_down_10kb (sum of reads in TSS_up_10kb & TES_down_10kb)
DensityTSS_TES=$( cat ${resultsFileRSeQC} | grep -iE "TSS_up_10kb|TES_down_10kb" | cut -f4 | awk '{SUM += $4} END {print SUM}' )

# write of these densities in the file ${samplename}_reads_distribution.txt
echo "Density in Exons (Tags/Kb): $DensityExons" >> ${readDistributionFile}
echo "Density in introns (Tags/Kb): $DensityIntrons" >> ${readDistributionFile}
echo "Density in TSS up 10kb / TES down 10kb (Tags/Kb): $DensityTSS_TES" >> ${readDistributionFile}

echo "****************************************"
echo "`date`: calculation of reads distribution : done"
echo "****************************************"
