#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# This script creates a table of all the statistics for each sample
# e.g sequencing quality statistics, tophat aligment statistics, reads distribution ...
# available in ${OUTDIR}/Quality/all_statistics_results.txt (statisticsFile)
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh OverallStatistics.sh -c config-file -s samplename -f fastq-file -r results-RSeQC -o OUTDIR [-h]

		Description:
		   This script creates a table of all the statistics for each sample

		Options:
			-c, --config-file Config file of the run.
		        This option is required. 
		        Generally named : config.sh
			-s, --samplename Samplename of the current file
		        This option is required.
			-f, --fastq-file a fastq file in .gz format, located in the Raw_data folder (or Processed_data if Spikes=yes)
		        This option is required.
			-r, --resultsRSeQC The results file from RSeQC, creates by the scripts ReadsDistribution.sh AND StrandSpecificity.sh
		        This option is required.
		 	-o, --outdir OUTDIR of the run
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
ARGS=$(getopt -o "c:s:f:r:o:h" --long "config-file:,samplename:,fastq-file:,resultsRSeQC:,outdir:,help" -- "$@" 2> /dev/null)
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
		-f|--fastq-file)
			FILE="$2"
			shift 2 
			;;
		-r|--resultsRSeQC)
			RESULTS_RSEQC="$2"
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
[ "$SAMPLENAME" == "" ] || [ "$CONFIG_FILE" == "" ] || [ "$FILE" == "" ] || [ "$RESULTS_RSEQC" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --samplename, --config-file, --fastq-file, --resultsRSeQC and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------

echo "****************************************************************"
echo "`date`: calculation of overall statistics : start"
echo "****************************************************************"

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

dir=${OUTDIR}/Alignment/${SAMPLENAME}
eval DESIGN_FILE=$DESIGN_FILE
eval STATISTICS_FILE=$STATISTICS_FILE
resultsFileRSeQC="${OUTDIR}/Quality/RSeQC/resultsTable_${SAMPLENAME}.txt"

# collecting values ​​of interest
sampleID=$( grep ${SAMPLENAME} ${DESIGN_FILE} | cut -f2 )
seqname=`zcat $FILE | head -1 | awk -F ":" '{print $1}' | awk -F "@" '{print $2}'`
flowcell=`zcat $FILE | head -1 | awk -F ":" '{print $3}'`
line=`zcat $FILE | head -1 | awk -F ":" '{print $4}'`
seqlengthplus1=`zcat $FILE | head -2 | tail -1 | wc | awk -F " " '{print $3}'`
seqlength=$(( ${seqlengthplus1} - 1 ))

if [ "$PAIRED" = "no" ]
then

	uniqueNumberSequences=$( grep "Unique number of sequences" $OUTDIR/Quality/stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	totalNumberSequences=$( grep "Total number of sequences" $OUTDIR/Quality/stats_$SAMPLENAME.txt | cut -d"|" -f 2 | sed "s/ //" )

	if [ "$SPIKES" = "yes" ]
	then

		totalNumberSequences_withoutSpikes=$(grep "Number of sequences without spikes" $OUTDIR/Quality/stats_$SAMPLENAME.txt | cut -d"|" -f 2 | sed "s/ //")
		numberSpikes=$( bc <<<"scale=0;(${totalNumberSequences}-${totalNumberSequences_withoutSpikes})" )
		percentSpikes=$( bc <<<"scale=3;(${numberSpikes}*100/${totalNumberSequences})" )

	fi

elif [ "$PAIRED" = "yes" ]
then

	# We have done a "paste" bewteen .R1 & .R2 to count unique sequences and we write this number only in the .R1 file to recover it
	uniqueNumberSequences=$( grep "Unique number of sequences" $OUTDIR/Quality/stats_$SAMPLENAME.R1.txt | cut -d"|" -f2 )
	# Averaging the total number of reads in both fastq files
	totalNumberSequences=$( grep "Total number of sequences" $OUTDIR/Quality/stats_$SAMPLENAME.R1.txt | cut -d"|" -f 2 | sed "s/ //" )

	if [ "$SPIKES" = "yes" ]
	then

		totalNumberSequences_withoutSpikes=$(grep "Number of sequences without spikes" $OUTDIR/Quality/stats_$SAMPLENAME.R1.txt | cut -d"|" -f 2 | sed "s/ //")
		numberSpikes=$( bc <<<"scale=0;(${totalNumberSequences}-${totalNumberSequences_withoutSpikes})" )
		percentSpikes=$( bc <<<"scale=3;(${numberSpikes}*100/${totalNumberSequences})" )

	fi

fi

total_unique=$( bc <<<"scale=3;(${totalNumberSequences}/${uniqueNumberSequences})" ) 

if [ "$SPIKES" = "yes" ]
then
	totalNumberSequencesToUse=${totalNumberSequences_withoutSpikes}
else
	totalNumberSequencesToUse=${totalNumberSequences}
fi

if [ "$PAIRED" = "yes" ]
then

	filtered=$( head -3 $dir/logs/prep_reads.log | tail -1 | cut -d" " -f1 )
	percent_filtered=$( bc <<<"scale=3;(${filtered}*100/${totalNumberSequencesToUse})" ) 

	alignedReads_by_Tophat2=$( grep "Nb of aligned pairs calculated by TopHat2" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentAlignedReads_by_Tophat2=$( bc <<<"scale=3;(${alignedReads_by_Tophat2}*100/${totalNumberSequencesToUse})" )

	uniqueAlignedReads_by_Tophat2=$( grep "Nb of uniquely aligned pairs calculated by TopHat2" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentUniqueAlignedReads_by_Tophat2=$( bc <<<"scale=3;(${uniqueAlignedReads_by_Tophat2}*100/${alignedReads_by_Tophat2})" )
	multipleAlignedReads_by_Tophat2=$( grep "Nb of multiple aligned pairs calculated by TopHat2" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentMultipleAlignedReads_by_Tophat2=$( bc <<<"scale=3;(${multipleAlignedReads_by_Tophat2}*100/${alignedReads_by_Tophat2})" )
	discordantPairsAlignment_by_Tophat2=$( grep "discordant pair alignments" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentdiscordantPairsAlignment_by_Tophat2=$( bc <<<"scale=3;(${discordantPairsAlignment_by_Tophat2}*100/${alignedReads_by_Tophat2})" )
	
	alignedReads=$( grep "Nb of aligned reads" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentAlignedReads=$( bc <<<"scale=3;(${alignedReads}*100/${totalNumberSequencesToUse})" )
	uniqueAlignedReads=$( grep "Nb of uniquely aligned reads" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentUniqueAlignedReads=$( bc <<<"scale=3;(${uniqueAlignedReads}*100/${alignedReads})" )
	multipleAlignedReads=$( grep "Nb of multiple aligned reads" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentMultipleAlignedReads=$( bc <<<"scale=3;(${multipleAlignedReads}*100/${alignedReads})" )

elif [ "$PAIRED" = "no" ]
then

	filtered=$( head -3 $dir/logs/prep_reads.log | tail -1 | cut -d" " -f1 )
	percent_filtered=$( bc <<<"scale=3;(${filtered}*100/${totalNumberSequencesToUse})" ) 

	alignedReads=$( grep "Nb of aligned reads" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentAlignedReads=$( bc <<<"scale=3;(${alignedReads}*100/${totalNumberSequencesToUse})" )

	uniqueAlignedReads=$( grep "Nb of uniquely aligned reads" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentUniqueAlignedReads=$( bc <<<"scale=3;(${uniqueAlignedReads}*100/${alignedReads})" )
	multipleAlignedReads=$( grep "Nb of multiple aligned reads" $OUTDIR/Alignment/tophat_stats_$SAMPLENAME.txt | cut -d"|" -f2 )
	percentMultipleAlignedReads=$( bc <<<"scale=3;(${multipleAlignedReads}*100/${alignedReads})" )

fi

# we do the sum of tag in exons (CDS+UTR), then we do the sum of total bases representing the exons and finaly we
# calculate the Density of exons (tags/kb)
TagExons=$( cat ${resultsFileRSeQC} | grep -iE "CDS|UTR" | cut -f3 | awk '{SUM += $3} END {print SUM}' )
TotalBasesExons=$( cat ${resultsFileRSeQC} | grep -iE "CDS|UTR" | cut -f2 | awk '{SUM += $2} END {print SUM}' )
DensityExons=$( bc <<<"scale=2;(${TagExons}*1000/${TotalBasesExons})" ) 
DensityIntrons=$( cat ${resultsFileRSeQC} | grep "Introns" | cut -f3 | awk '{SUM += $4} END {print SUM}' )
TagTSS_TES=$( cat ${resultsFileRSeQC} | grep -iE "TSS_up_10kb|TES_down_10kb" | cut -f3 | awk '{SUM += $3} END {print SUM}' )
TotalBasesTSS_TES=$( cat ${resultsFileRSeQC} | grep -iE "TSS_up_10kb|TES_down_10kb" | cut -f2 | awk '{SUM += $2} END {print SUM}' )
DensityTSS_TES=$( bc <<<"scale=2;(${TagTSS_TES}*1000/${TotalBasesTSS_TES})" ) 
assigned=`grep '\(ERCC\|ENS\)' $OUTDIR/Alignment/${SAMPLENAME}/${SAMPLENAME}_htseq.txt | cut -f 2 | perl -e '$sum=0;while ($line=<>) {$sum=$sum+$line;}; print $sum; '`
percentAssigned=$( bc <<<"scale=3;(${assigned}*100/${uniqueAlignedReads})" ) 
nofeature=`grep no_feature $OUTDIR/Alignment/${SAMPLENAME}/${SAMPLENAME}_htseq.txt | cut -f 2`
percentNofeature=$( bc <<<"scale=3;(${nofeature}*100/${uniqueAlignedReads})" ) 
ambiguous=`grep ambiguous $OUTDIR/Alignment/${SAMPLENAME}/${SAMPLENAME}_htseq.txt | cut -f 2`
percentAmbiguous=$( bc <<<"scale=3;(${ambiguous}*100/${uniqueAlignedReads})" ) 
notunique=`grep alignment_not_unique $OUTDIR/Alignment/${SAMPLENAME}/${SAMPLENAME}_htseq.txt | cut -f 2`

if [ "$STRANDED" = "reverse" ]
then

	if [ "$PAIRED" = "no" ]
	then

		percentNumberForward=$( grep "++,--" ${OUTDIR}/Quality/RSeQC/strand_specificity/${SAMPLENAME}_strand_specificity.txt | cut -d: -f2 )
		percentNumberReverse=$( grep "+-,-+" ${OUTDIR}/Quality/RSeQC/strand_specificity/${SAMPLENAME}_strand_specificity.txt | cut -d: -f2 )
	
	elif [ "$PAIRED" = "yes" ]
	then

		percentNumberForward=$( grep "1++,1--,2+-,2-+" ${OUTDIR}/Quality/RSeQC/strand_specificity/${SAMPLENAME}_strand_specificity.txt | cut -d: -f2 )
		percentNumberReverse=$( grep "1+-,1-+,2++,2--" ${OUTDIR}/Quality/RSeQC/strand_specificity/${SAMPLENAME}_strand_specificity.txt | cut -d: -f2 )

	else

		print_err_and_exit "What is your sequencing protocol ? Paired-end = yes or no ?"

	fi

fi



## write of the statistics in the ${STATISTICS_FILE} file.

echo -ne "\n${SAMPLENAME}\t${sampleID}\t${seqname}\t${flowcell}\t${line}\t${seqlength}\t" >> ${STATISTICS_FILE}

if [ "$SPIKES" = "yes" ]
then

	echo -ne "${totalNumberSequences}\t${totalNumberSequences_withoutSpikes}\t${uniqueNumberSequences}\t${total_unique}\t${percentSpikes}\t" >> ${STATISTICS_FILE}

else 

	echo -ne "${totalNumberSequences}\t${uniqueNumberSequences}\t${total_unique}\t" >> ${STATISTICS_FILE}

fi

if [ "$PAIRED" = "no" ]
then

	echo -ne "${filtered}\t${percent_filtered}\t${alignedReads}\t${percentAlignedReads}\t${uniqueAlignedReads}\t${percentUniqueAlignedReads}\t${multipleAlignedReads}\t${percentMultipleAlignedReads}\t" >> ${STATISTICS_FILE}

else

	echo -ne "${filtered}\t${percent_filtered}\t${alignedReads_by_Tophat2}\t${percentAlignedReads_by_Tophat2}\t${uniqueAlignedReads_by_Tophat2}\t${percentUniqueAlignedReads_by_Tophat2}\t${multipleAlignedReads_by_Tophat2}\t${percentMultipleAlignedReads_by_Tophat2}\t${discordantPairsAlignment_by_Tophat2}\t${percentdiscordantPairsAlignment_by_Tophat2}\t${alignedReads}\t${percentAlignedReads}\t${uniqueAlignedReads}\t${percentUniqueAlignedReads}\t${multipleAlignedReads}\t${percentMultipleAlignedReads}\t" >> ${STATISTICS_FILE}

fi

if [ "$STRANDED" = "reverse" ]
then

	echo -ne "${DensityExons}\t${DensityIntrons}\t${DensityTSS_TES}\t${assigned}\t${percentAssigned}\t${nofeature}\t${percentNofeature}\t${ambiguous}\t${percentAmbiguous}\t${notunique}\t${percentNumberForward}\t${percentNumberReverse}" >> ${STATISTICS_FILE}

else

	echo -ne "${DensityExons}\t${DensityIntrons}\t${DensityTSS_TES}\t${assigned}\t${percentAssigned}\t${nofeature}\t${percentNofeature}\t${ambiguous}\t${percentAmbiguous}\t${notunique}" >> ${STATISTICS_FILE}

fi

echo "****************************************"
echo "`date`: calculation of overall statistics : done"
echo "****************************************"
