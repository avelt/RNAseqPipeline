#! /bin/sh
#--------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#--------------------------------------------------------------------------------------------------------------------------------------------
# This script is dedicated to the separation of the spike's reads and the other reads
# together with the processing of the RNAseq sample (ie all the reads without spikes).
# Spikes will be processed separately by a dedicated script.
#--------------------------------------------------------------------------------------------------------------------------------------------
# if spikes = yes, this script is executed. It separates the raw fastq file into two files:
# - one fastq file without reads corresponding to spikes (and we have called "${samplename}_without_Spikes.fastq" in run.sh)
# - one sam file containing all reads of spikes
# Then it execute the scripts "QualityAnalysis.sh" and "RNAseqPipeline1sample.sh" on the fastq file which doesn't contain the spikes' reads.
# (whereas if spikes = no, these two scripts are executed directly on the raw fastq file through the script "run.sh").
#--------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh AnalysisAllSampleWithoutSpikes.sh -rf raw-fastq-file -pf processed-fastq-file -rn raw-name -rs raw-samplename -c config_file -o outdir [-h]

		Description:
		    This script is dedicated to the separation of the spike's reads and the other reads together with the processing of the RNAseq 
		    sample (ie all the reads without spikes). Spikes will be processed separately by a dedicated script.

		Options:
		 	-rf, --raw-fastq raw fastq file in .gz format, located in the rawdata folder
		        This option is required. 
		        Generally of type : /path/to/Num_project_Name_author/rawdata
		 	-pf, --processed-fastq fastq file in .gz format, containing reads which are not from Spikes, and located in the Processed_data folder
		        This option is required. 
		        Generally of type : /path/to/Num_project_Name_author/Date/Processed_data
		 	-rn, --raw-name Name of the raw fastq file, defined in run.sh
		        This option is required. 
		        Generally of type : /path/to/name_fastq_file
		 	-rs, --raw-samplename Samplename of the raw fastq file, defined in run.sh
		        This option is required. 
		        Generally of type : name_fastq_file
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

#------------------------------------------------------------------------------------------------------------------
# Getting parameters from the input
#------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "rf:pf:rn:rs:c:o:h" --long "raw-fastq:,processed-fastq:,raw-name:,raw-samplename:,config-file:,outdir:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-rf|--raw-fastq)
			RAW_FASTQ="$2"
			shift 2 
			;;
		-pf|--processed-fastq)
			PROCESSED_FASTQ="$2"
			shift 2 
			;;
		-rn|--raw-name)
			RAW_NAME="$2"
			shift 2 
			;;
		-rs|--raw-samplename)
			RAW_SAMPLENAME="$2"
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
[ "$RAW_FASTQ" == "" ] || [ "$PROCESSED_FASTQ" == "" ] || [ "$RAW_NAME" == "" ] || [ "$RAW_SAMPLENAME" == "" ] || [ "$CONFIG_FILE" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --raw-fastq, --processed-fastq, --raw-name, --raw-samplename, --config-file and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------

# to wait until the script ends.
touch $OUTDIR/wait_run_end/${RAW_SAMPLENAME}_tmp.txt

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

#path/to/file without spike's reads
file_without_spikes_in_GZ_format=${PROCESSED_FASTQ}
file_without_spikes=$( echo $file_without_spikes_in_GZ_format | sed "s/\.gz//" )

#path/to/sample_name (= file name without final .fastq.gz)
# there is the command "sed -e 's/\.R1//'" because when we work on paired-end data, we send only the .R1 fastq file to this script, in order to recover the samplename.
name_file_without_spikes=$( get_name ${file_without_spikes_in_GZ_format} )

#sample name
samplename_file_without_spikes=$( get_sample_name ${file_without_spikes_in_GZ_format} )

# alignment of spikes and separation in two FASTQ files to perform all 
# the analysis of spikes separately from other reads.
make_directories ${OUTDIR}/Alignment/${samplename_file_without_spikes}	
#gunzip ${RAW_NAME}*fastq.gz
wait

#------------------------------------------------------------------------------------------------------------------
# for this alignment step :
# we use another sbatch command because Bowtie uses a number of cpu, so send this command on another sbatch allows 
# to cleanly assign the job on the different cluster.
#------------------------------------------------------------------------------------------------------------------

echo " `date`: alignment of spikes and separation in two FASTQ files "

if [ "$PAIRED" = "no" ]
then

	touch ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam 
	mkfifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo
	zcat ${RAW_FASTQ} > ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo &

	#raw_fastq_gunzip=$( echo ${RAW_FASTQ} | sed 's/.gz//' )
	echo "${RAW_SAMPLENAME}" >> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt

	if echo ${BOWTIE_VERSION} | grep -qE "bowtie-1|bowtie-0"
	then

		# with bowtie1, the option "--un" allows to recover the unaligned reads in a separate fastq file
		# corresponding to reads which are not from Spikes.
		# the .sam file will contain the bowtie alignment, ie all the reads corresponding to spikes
			
		echo "#------------------------------------------"
		echo "Bowtie commands : Alignment on Spikes =>"
		echo "${BIN}/${BOWTIE_VERSION}/bowtie -p ${NBPROC} --un ${file_without_spikes} ${INDEX_SPIKES} ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt"
		echo "#------------------------------------------"

		${BIN}/${BOWTIE_VERSION}/bowtie -p ${NBPROC} --un ${file_without_spikes} ${INDEX_SPIKES} ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt

	elif echo ${BOWTIE_VERSION} | grep -q "bowtie-2"
	then
				
		# with bowtie 2, outputs of the alignment are different of bowtie1.
		# the option "--un" allows to recover the unaligned reads in a separate fastq file
		# BUT, the output .sam file contains all the reads, those aligned as well as unaligned
		# We remove the unaligned reads of the .sam file, ie we remove the lines that don't contain the "AS" flag
		# in order to kepp the spikes only.
		# (there is a "--no-unal" option for that the unaligned reads don't appear in the .sam file 
		# BUT the "--un-conc" option no longer works and fastq file is empty)
		
		echo "#------------------------------------------"
		echo "Bowtie command : Alignment on Spikes =>"
		echo "${BIN}/${BOWTIE_VERSION}/bowtie2 -U ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo -p ${NBPROC} --un ${file_without_spikes} ${INDEX_SPIKES} -S ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt &"
		echo "#------------------------------------------"

		${BIN}/${BOWTIE_VERSION}/bowtie2 -U ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo -p ${NBPROC} --un ${file_without_spikes} ${INDEX_SPIKES} -S ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt &

	else

		print_err_and_exit "Bowtie version is not known. Exit RNAseqPipeline."

	fi

	wait
	rm ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}.fifo
	gzip ${file_without_spikes}
	# removal of reads do not correspond to spikes
	sed -i '/AS:*/!d' ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam
	wait

	# QualityAnalysis.sh		
	
	echo "#------------------------------------------"
	echo "Quality analysis command :"
	echo "bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${file_without_spikes_in_GZ_format} --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &"
	echo "#------------------------------------------"

	bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${file_without_spikes_in_GZ_format} --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &

	# RNAseqPipeline1sample.sh
	
	echo "#------------------------------------------"
	echo "RNAseqPipeline1sample on $RAW_SAMPLENAME :"
	echo "sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh --fastq-file ${file_without_spikes_in_GZ_format} --config-file ${CONFIG_FILE} --outdir $OUTDIR >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &"
	echo "#------------------------------------------"

	sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh --fastq-file ${file_without_spikes_in_GZ_format} --config-file ${CONFIG_FILE} --outdir $OUTDIR >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &

	wait
	
elif [ "$PAIRED" = "yes" ]
then

	touch ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam 
	wait
	mkfifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo
	mkfifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo
	zcat ${RAW_NAME}.R1.fastq.gz > ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo &
	zcat ${RAW_NAME}.R2.fastq.gz > ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo &

	if echo ${BOWTIE_VERSION} | grep -qE "bowtie-1|bowtie-0"
	then
					
		# with bowtie1, the option "--un" allows to recover the unaligned reads in a separate fastq file
		# corresponding to reads which are not from Spikes.
		# the .sam file will contain the bowtie alignment, ie all the reads corresponding to spikes
		echo "${RAW_SAMPLENAME}" >> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt
		
		echo "#------------------------------------------"
		echo "Bowtie commands : Alignment on Spikes =>"
		echo "${BIN}/${BOWTIE_VERSION}/bowtie -p ${NBPROC} --un ${name_file_without_spikes}.fastq ${INDEX_SPIKES} -1 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo -2 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt &"

		${BIN}/${BOWTIE_VERSION}/bowtie -p ${NBPROC} --un ${name_file_without_spikes}.fastq ${INDEX_SPIKES} -1 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo -2 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt &

	elif echo ${BOWTIE_VERSION} | grep -q "bowtie-2"
	then
				
		# with bowtie 2, outputs of the alignment are different of bowtie1.
		# the option "--un" allows to recover the unaligned reads in a separate fastq file
		# BUT, the output .sam file contains all the reads, those aligned as well as unaligned
		# We remove the unaligned reads of the .sam file, ie we remove the lines that don't contain the "AS" flag
		# in order to kepp the spikes only.
		# (there is a "--no-unal" option for that the unaligned reads don't appear in the .sam file 
		# BUT the "--un-conc" option no longer works and fastq file is empty)
		
		echo "${RAW_SAMPLENAME}" >> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt

		echo "#------------------------------------------"
		echo "Bowtie commands : Alignment on Spikes =>"
		echo "${BIN}/${BOWTIE_VERSION}/bowtie2 -1 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo -2 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo -p ${NBPROC} --un-conc ${name_file_without_spikes}.R%.fastq ${INDEX_SPIKES} -S ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt &"

		${BIN}/${BOWTIE_VERSION}/bowtie2 -1 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo -2 ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo -p ${NBPROC} --un-conc ${name_file_without_spikes}.R%.fastq ${INDEX_SPIKES} -S ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam &>> $OUTDIR/Quality/Spikes/${RAW_SAMPLENAME}_stats_Spikes.txt &
			
	else

		print_err_and_exit "Bowtie version is not known. Exit RNAseqPipeline."

	fi

	wait
	
	# the output of bowtie1 is ${name}_without_Spikes.1.fastq and the output of bowtie2 is
	# ${name}_without_Spikes-1.fastq, so we rename to have the same filename regardless of the bowtie version.
	if echo ${BOWTIE_VERSION} | grep -qE "bowtie-1|bowtie-0"
	then
		if [ ! -f ${name_file_without_spikes}.R1.fastq ]
		then

			mv ${name_file_without_spikes}.1.fastq ${name_file_without_spikes}.R1.fastq || echo "Cannot move ${name_file_without_spikes}*1.fastq to ${name_file_without_spikes}.R1.fastq  - EXIT " && exit 1;
			mv ${name_file_without_spikes}.2.fastq ${name_file_without_spikes}.R2.fastq || echo "Cannot move ${name_file_without_spikes}*2.fastq to ${name_file_without_spikes}.R2.fastq  - EXIT " && exit 1;

		fi
	fi

	wait

	rm ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R1.fifo ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_R2.fifo
	gzip ${name_file_without_spikes}.R*.fastq

	sed -i '/AS:*/!d' ${OUTDIR}/Alignment/${samplename_file_without_spikes}/${RAW_SAMPLENAME}_Spikes.sam
	
	wait

	# two quality analyses are performed, one for .R1 file and one for .R2 file.

	# RNAseqPipeline1sample.sh
	echo "#------------------------------------------"
	echo "RNAseqPipeline1sample on $samplename_file_without_spikes :"
	echo "sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh --fastq-file ${name_file_without_spikes}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir $OUTDIR >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &"
	echo "#------------------------------------------"

	sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh --fastq-file ${name_file_without_spikes}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir $OUTDIR >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &

	# QualityAnalysis.sh on .R1 file	
	echo "#------------------------------------------"
	echo "Quality analysis command :"
	echo "bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name_file_without_spikes}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes}"
	echo "#------------------------------------------"

	bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name_file_without_spikes}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes}

	wait

	# QualityAnalysis.sh on .R2 file		
	echo "#------------------------------------------"
	echo "Quality analysis command :"
	echo "bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name_file_without_spikes}.R2.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &"
	echo "#------------------------------------------"

	bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name_file_without_spikes}.R2.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename_file_without_spikes} &

	wait

else

	print_err_and_exit "The paired-end variable is misnamed. Exit RNAseqPipeline."

fi	
		
# to say that the script ends
rm $OUTDIR/wait_run_end/${RAW_SAMPLENAME}_tmp.txt
