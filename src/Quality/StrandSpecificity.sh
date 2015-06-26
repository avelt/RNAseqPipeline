#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# Analysis of strand specificity in the experiment
# The script "strand_specificity.py" generates a file with
# the number of reads on each strand or others features (intergenic, overlap, unknown ...)
# Requires the existence of the gtf reference file in bed format
# If this bed file does not exist, it is automatically created with the script gtf2bed.pl
# Requires fasta file of the genome use. For this, we use $GENOME_VERSION value to convert in corresponding fasta file
# the conversion file is in RNAseqReport/Files folder, and is called "genome_conversion.txt"
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh StrandSpecificity.sh -c config-file -s samplename -o outdir [-h]

		Description:
		    Analysis of strand specificity in the experiment.

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
echo "`date`: calculation of strand specificity : start"
echo "****************************************************************"

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

cd ${BIN}/${RSEQC_VERSION}
export PYTHONPATH=${PYTHONPATH}/${RSEQC_PYTHONPATH}/

resultsFileRSeQC=$OUTDIR/Quality/RSeQC/resultsTable_${SAMPLENAME}.txt

bedfile=$( create_bed_from_gtf ${GTF_FILE} ${SCRIPTSDIR} )

# creation of file containing strand specificity information
strandSpecificityFile=${OUTDIR}/Quality/RSeQC/strand_specificity/${SAMPLENAME}_strand_specificity.txt
touch $strandSpecificityFile

# if there are no spikes, stranded specificity is tested on all reads
if [ "$SPIKES" = "no" ]
then
	
		# recovery of genome version use (hg19, mm9, hg19spikes, mm9spikes ...)
		genomeAssembly=$( echo ${GENOME_VERSION} | grep -o "[^/]*$" )
		# recovery of fasta file corresponding to genome version use
		FastaFile=$( grep ${genomeAssembly} ${REPORTDIR}/Files/genome_conversion.txt | cut -f3 )
		
		echo "python -u infer_experiment.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -r ${bedfile} &>> ${resultsFileRSeQC}"
		
		# if we use furious, we use softwareSL and python 2.7 must be recovered (else python 2.6 is used bu default)
		echo $BIN | if grep -q "softwareSL\|biopuces"; then
		
			${BIN}/Python/bin/python -u infer_experiment.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -r ${bedfile} &>> ${resultsFileRSeQC}
		
		else
		
			python -u infer_experiment.py -i ${OUTDIR}/Alignment/${SAMPLENAME}/${SAMPLENAME}_alignment.bam -r ${bedfile} &>> ${resultsFileRSeQC}
		
		fi
		
		if [ "$PAIRED" = "no" ]
		then

			#reads count on forward strand
			countForwardStrand=$( grep 'Fraction of reads explained by "++,--":' ${resultsFileRSeQC} | cut -d: -f2 )
			#reads count on reverse strand
			countReverseStrand=$( grep 'Fraction of reads explained by "+-,-+":' ${resultsFileRSeQC} | cut -d: -f2 )
			#other reads
			countOther=$( grep 'Fraction of reads failed to determine:' ${resultsFileRSeQC} | cut -d: -f2 )
			
			# calculation of percentage and write to file ${samplename}_strand_specificity.txt
			echo "Fraction of reads explained by '++,--': $countForwardStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads explained by '+-,-+': $countReverseStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads failed to determine: $countOther" >> ${strandSpecificityFile}
		
		elif [ "$PAIRED" = "yes" ]
		then

			#reads count on forward strand
			countForwardStrand=$( grep 'Fraction of reads explained by "1++,1--,2+-,2-+":' ${resultsFileRSeQC} | cut -d: -f2 )
			#reads count on reverse strand
			countReverseStrand=$( grep 'Fraction of reads explained by "1+-,1-+,2++,2--":' ${resultsFileRSeQC} | cut -d: -f2 )
			#other reads
			countOther=$( grep 'Fraction of reads failed to determine:' ${resultsFileRSeQC} | cut -d: -f2 )
			
			# calculation of percentage and write to file ${samplename}_strand_specificity.txt
			echo "Fraction of reads explained by '1++,1--,2+-,2-+' (forward) : $countForwardStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads explained by '1+-,1-+,2++,2--' (reverse): $countReverseStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads failed to determine: $countOther" >> ${strandSpecificityFile}

		else

			print_err_and_exit "What is your sequencing protocol ? Paired-end = yes or no ? Exit RNAseqPipeline"

		fi
		
else 

		samFileSpikes=`ls ${OUTDIR}/Alignment/${SAMPLENAME}/*_Spikes.sam`
		
		if grep -q "@SQ" ${samFileSpikes}
		then
		
			echo " ${samFileSpikes} header OK "
			
		else
		
			dirName=$( dirname ${samFileSpikes} )
			sampleName=$( basename ${samFileSpikes} | sed "s/Spikes/Spikes_for_strandSpecificity/" )
			samFileSpikesOK="${dirName}/${sampleName}"
			cat ${SCRIPTSDIR_SPIKES}/Files/SAM_header_spikes.txt ${samFileSpikes} > ${samFileSpikesOK}
			samFileSpikes=${samFileSpikesOK}
		
		fi

		echo "python -u infer_experiment.py -i ${samFileSpikes} -r ${SPIKES_ANNOTATIONS} &>> ${resultsFileRSeQC}"
		
		# if we use furious, we use softwareSL and python 2.7 must be recovered (else python 2.6 is used bu default)
		echo $BIN | if grep -q "softwareSL\|biopuces"; then
		
			${BIN}/Python/bin/python -u infer_experiment.py -i ${samFileSpikes} -r ${SPIKES_ANNOTATIONS} &>> ${resultsFileRSeQC}
	
		else
		
			python -u infer_experiment.py -i ${samFileSpikes} -r ${SPIKES_ANNOTATIONS} &>> ${resultsFileRSeQC}
		
		fi

		if [ "$PAIRED" = "no" ]
		then
	
			#reads count on forward strand
			countForwardStrand=$( grep '++,--' ${resultsFileRSeQC} | cut -d: -f2 )
			#reads count on reverse strand
			countReverseStrand=$( grep '+-,-+' ${resultsFileRSeQC} | cut -d: -f2 )
			#other reads
			countOther=$( grep 'Fraction of reads failed to determine:' ${resultsFileRSeQC} | cut -d: -f2 )
			
			# calculation of percentage and write to file ${samplename}_strand_specificity.txt
			echo "Fraction of reads explained by '++,--': $countForwardStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads explained by '+-,-+': $countReverseStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads failed to determine: $countOther" >> ${strandSpecificityFile}

		elif [ "$PAIRED" = "yes" ]
		then

			#reads count on forward strand
			countForwardStrand=$( grep '1++,1--,2+-,2-+' ${resultsFileRSeQC} | cut -d: -f2 )
			#reads count on reverse strand
			countReverseStrand=$( grep '1+-,1-+,2++,2--' ${resultsFileRSeQC} | cut -d: -f2 )
			#other reads
			countOther=$( grep 'Fraction of reads failed to determine:' ${resultsFileRSeQC} | cut -d: -f2 )
			
			# calculation of percentage and write to file ${samplename}_strand_specificity.txt
			echo "Fraction of reads explained by '1++,1--,2+-,2-+' (forward): $countForwardStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads explained by '1+-,1-+,2++,2--' (reverse): $countReverseStrand" >> ${strandSpecificityFile}
			echo "Fraction of reads failed to determine: $countOther" >> ${strandSpecificityFile}

		else

			print_err_and_exit "What is your sequencing protocol ? Paired-end = yes or no ? Exit RNAseqPipeline"

		fi
		
fi
	
echo "`date`: calculation of strand specificity : done"

	
