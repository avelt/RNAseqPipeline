#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Céline Keime / Amandine VELT
# Mail: keime@igbmc.fr / velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
# This script is the main script of the RNA-seq pipeline. Use it to launch the whole pipeline.
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
#
# PREPARATION BEFORE EXECUTING THE RNA-SEQ PIPELINE
#
###################################################################################################################


###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    nohup sh run.sh -o outdir [-h] > \$outdir/nohup.allsamplesanalysis

		Description:
		    run.sh is the main script of the RNA-seq analysis pipeline of the IGBMC Microarray and Sequencing platform.
		    Don't forget the following :
		    - fill in correctly the config.sh file and place it in \$outdir before running run.sh script.
		    - create a path/to/outdir/rawdata folder containing all the raw fastq files in .gz format.
		    - create a \$outdir/design_file.txt file (see documentation to know how write it)

		Options:
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
ARGS=$(getopt -o "o:h" --long "outdir:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-o|--outdir)
			OUTDIR_RUN="$2"
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
[ ! -d "$OUTDIR_RUN" ] && \
	echo "$OUTDIR_RUN directory doesn't exist ! Rerun with valid directory path." && exit 1;
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Initialization and checking
#------------------------------------------------------------------------------------------------------------------

# We check that config file exists
# this file contains all the parameters of the run and must be filled before running the RNAseqPipeline.

if [ ! -f "$OUTDIR_RUN/config.sh" ]
then
	echo "$OUTDIR_RUN/config.sh file doesn't exist. Create it and re-run RNA-seq pipeline." 1>&2
	exit 1 
fi

if ! source "$OUTDIR_RUN/config.sh"
then
	echo "Cannot load $OUTDIR_RUN/config.sh, create it and re-run RNA-seq pipeline." 1>&2
	exit 1 
else 
	# Getting parameters
	source "$OUTDIR_RUN/config.sh"
fi

# We check that utils file exists
# this file contains all the functions used by the RNAseqPipeline and must be accessible
if ! source "$SCRIPTSDIR/src/Utilities/utils.sh"
then
	echo "$SCRIPTSDIR/src/Utilities/utils.sh don't find. The $SCRIPTSDIR path is ok? Check and re-run."
	exit 1
else
	source $SCRIPTSDIR/src/Utilities/utils.sh
fi

# We create one folder by RNAseqPipeline's run. 
# So, here, we create, in $Outdir, a folder named with the date fo the run.
OUTDIR=$( echo ${OUTDIR_RUN}/${OUTDIR_ANALYSIS} )

# this function create the specific outdir + the directory which will contain all the scripts for this specific outdir
make_directories $OUTDIR $OUTDIR/Scripts

cp $OUTDIR_RUN/config.sh $OUTDIR/Scripts/config.sh
CONFIG_FILE="$OUTDIR/Scripts/config.sh"

# we check if the rawdata folder exists
# otherwise, we exit the pipeline
if [ ! -d "$OUTDIR_RUN/${DATA_FOLDER}" ]
then

	print_err_and_exit "Cannot find $OUTDIR_RUN/${DATA_FOLDER} which contain all the raw fastq files. Create it and re-run RNA-seq pipeline." 1>&2

else

    DATADIR="$OUTDIR_RUN/${DATA_FOLDER}"

fi

# we copy all the RNAseqPipeline's scripts in the $OUTDIR/Scripts/ folder, in order to keep track of used scripts 
# since we are lead to make changes in the RNAseqPipeline
if [ ! -d "$OUTDIR/Scripts/src" ]
then

	cp -R $SCRIPTSDIR/src $OUTDIR/Scripts/src
	cp -R $SCRIPTSDIR/run.sh $OUTDIR/Scripts/run.sh

else 

	print_err_and_exit "$OUTDIR/Scripts/src directory exist, don't copy it from $SCRIPTSDIR/src. RNA-seq Pipeline continue." 1>&2

fi

# when you run the RNAseqPipeline, it is necessary to have the design_file corresponding to your samples
# if STATISTICAL_ANALYSIS='no', you just need the concordance between Sample_ID and Sample_name.
# if STATISTICAL_ANALYSIS='yes', you must write the design of the experiment (replicates, contrast, ...)

if [ ! -f "$OUTDIR_RUN/design_file.txt" ]
then

	print_err_and_exit "Cannot find $OUTDIR_RUN/design_file.txt which contain the design of the experiment. Create it and re-run RNA-seq pipeline. See the RNAseqPipeline's documentation to see how write this file." 1>&2
	
else

	cp "$OUTDIR_RUN/design_file.txt" "$OUTDIR/Scripts/design_file.txt"
	DESIGN_FILE="$OUTDIR/Scripts/design_file.txt"
	DESIGN_FILE_REPORT="$OUTDIR/Scripts/design_file_report.txt"
	cp $DESIGN_FILE $DESIGN_FILE_REPORT

fi

#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Attribution of good values to some variables of config file:
# because in config file, these variables have the 'text' $OUTDIR
# here, they take the 'value' of $OUTDIR, not define in config file !
#------------------------------------------------------------------------------------------------------------------

eval STATISTICS_FILE=$STATISTICS_FILE
eval GENE_BIOTYPE_FILE=$GENE_BIOTYPE_FILE
eval LOGFILE_SPIKES=$LOGFILE_SPIKES
eval LOGFILE_REPORT=$LOGFILE_REPORT
eval LOGFILE_WIG_TDF=$LOGFILE_WIG_TDF

#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
#
#This script launches RNAseqPipeline1sample.sh for all fastq.gz files located in $DATADIR
#
#------------------------------------------------------------------------------------------------------------------
cd $OUTDIR  || echo "Cannot go to $OUTDIR - Folder doesn't exist ?"
# for verification
echo "##############
Slurm variables :
"
echo "USE_SLURM : $USE_SLURM"
echo "SLURM_PARTITION : $SLURM_PARTITION"

echo "##############
Information relative to the project :
"
echo "DATADIR : $DATADIR"
echo "OUTDIR : $OUTDIR"
echo "PROJECT_NUMBER : $PROJECT_NUMBER"
echo "DESIGN_FILE : $DESIGN_FILE"

echo "##############
Information on library preparation and sequencing :
"
echo "STRANDED : $STRANDED"
echo "SPIKES : $SPIKES"
echo "PAIRED-END : $PAIRED"

echo "##############
Information on species, genome and annotations :
"
echo "SPECIES : $SPECIES"
echo "GENOME_VERSION : $GENOME_VERSION"
echo "GENOME_PATH : $GENOME_PATH"
echo "SPIKES_PATH : $SPIKES_PATH"
echo "BOWTIE_INDEXES : $BOWTIE_INDEXES"
echo "GTF_FILE : $GTF_FILE"

echo "##############
Options for the analyses performed by the pipeline :
"
echo "NBPROC : $NBPROC"
echo "STATISTICAL_ANALYSIS : $STATISTICAL_ANALYSIS"
echo "DIST_METHOD : $DIST_METHOD"
echo "HCLUST_METHOD : $HCLUST_METHOD"
echo "THRESHOLD_LOG_FC : $THRESHOLD_LOG_FC"
echo "THRESHOLD_ADJ_PVAL : $THRESHOLD_ADJ_PVAL"
echo "FIT_TYPE : $FIT_TYPE"
echo "WINDOW_SIZE : $WINDOW_SIZE"

echo "##############
Scripts and software versions and location:
"
echo "BIN : $BIN"
echo "RBIN : $RBIN"
echo "SCRIPTSDIR : $SCRIPTSDIR"
echo "UTILITIESDIR : $UTILITIESDIR"
echo "FASTQC_VERSION : $FASTQC_VERSION"
echo "TOPHAT_VERSION : $TOPHAT_VERSION"
echo "BOWTIE_VERSION : $BOWTIE_VERSION"
echo "SAMTOOLS_VERSION : $SAMTOOLS_VERSION"
echo "HTSEQ_VERSION : $HTSEQ_VERSION"
echo "RSEQC_VERSION : $RSEQC_VERSION"
echo "IGV_VERSION : $IGV_VERSION"
echo "PYTHONPATH : $PYTHONPATH"

echo "##############
Log files location:
"
echo "LOGS directory : `dirname $LOGFILE_REPORT`"
echo
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
#Create directories to store all results
#------------------------------------------------------------------------------------------------------------------

# creation of different outputs :
#output directory for all quality analysis files
#output directory for fastqc result files
#output directory for alignment result files
#output directory for results generated by RSeQC
#output directory for reads distribution results files
#output directory for gene body coverage
# output directory for log files

make_directories "$OUTDIR/Quality" \
                 "$OUTDIR/Quality/fastqc" \
                 "$OUTDIR/Alignment" \
                 "$OUTDIR/Quality/RSeQC" \
                 "$OUTDIR/Quality/RSeQC/reads_distribution" \
                 "$OUTDIR/Quality/RSeQC/geneBodyCoverage" \
                 "$OUTDIR/Quantification" \
                 "$OUTDIR/Logs"

if [ "$PAIRED" = "yes" ]
then
    make_directories "$OUTDIR/Quality/RSeQC/innerDistance" #output directory for inner distance calculation
fi

if [ "$SPIKES" = "yes" ]
then
    #output directory for spikes results files
    make_directories "$OUTDIR/Quality/Spikes/" \
                     "$OUTDIR/Processed_data/"
fi

if [ "$STRANDED" = "reverse" ]
then
    # only for STRANDED protocols
    make_directories "$OUTDIR/Quality/RSeQC/strand_specificity" #output directory for strand specificity results files
fi

# directory where all the temporary files are stored, 
# in order to wait that all the sbatch jobs are finished before continuing RNAseqPipeline.
make_directories $OUTDIR/wait_run_end
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# OVERALL STATISTICS FOR EACH SAMPLE
# create a table of all the statistics for each sample
# e.g sequencing quality statistics, tophat aligment statistics, reads distribution ...
# available in ${OUTDIR}/quality_analysis/all_statistics_results.txt (statisticsFile)
# This table is filled via RNAseqPipeline1sample.sh script
#------------------------------------------------------------------------------------------------------------------

touch ${STATISTICS_FILE}

# header of overall statistics table
echo -ne "Sample ID\tSample name\tSequencer\tFlow cell\tLine\tRead Length\t" >> ${STATISTICS_FILE}

# there are some supplementary colums depending on the options
# if spikes=yes, there are some colums on spike statistics
if [ "$SPIKES" = "yes" ]
then
    echo -ne "Total number of sequences (pairs for PE data)\tTotal number of sequences without spikes (pairs for PE data)\tUnique number of sequences (pairs for PE data)\tTotal/Unique number of sequences\t% of Spikes\t" >> ${STATISTICS_FILE}
else 
    echo -ne "Total number of sequences (pairs for PE data)\tUnique number of sequences (pairs for PE data)\tTotal/Unique number of sequences\t" >> ${STATISTICS_FILE}
fi

# if paired=yes, there are some supplementary columns (eg discordant pairs)
if [ "$PAIRED" = "no" ]
then
    echo -ne "Number of reads filtered out by Tophat\t% of reads filtered out by Tophat\tNumber of aligned reads\t% of aligned reads\tNumber of uniquely aligned reads\t% of uniquely aligned reads\tNumber of multiple aligned reads\t% of multiple aligned reads\t" >> ${STATISTICS_FILE}
else
    echo -ne "Number of reads filtered out by Tophat\t% of reads filtered out by Tophat\tNumber of aligned reads calculated by TopHat2\t% of aligned reads calculated by TopHat2\tNumber of uniquely aligned reads calculated by TopHat2\t% of uniquely aligned reads calculated by TopHat2\tNumber of multiple aligned reads calculated by TopHat2\t% of multiple aligned reads calculated by TopHat2\tNumber of discordant pairs alignment calculated by TopHat2\t% of discordant pairs alignment calculated by TopHat2\tNumber of aligned reads\t% of aligned reads\tNumber of uniquely aligned reads\t% of uniquely aligned reads\tNumber of multiple aligned reads\t% of multiple aligned reads\t" >> ${STATISTICS_FILE}
fi

# if stranded=reverse, there are two supplementary columns to see if the data are well reverse stranded.
if [ "$STRANDED" = "reverse" ]
then
    echo -ne "Density in exons (Tags/Kb)\tDensity in introns (Tags/Kb)\tDensity in TSS up 10kb / TES down 10kb (Tags/Kb)\tNumber of assigned reads\t% of assigned reads\tNumber of no feature reads\t% of no feature reads\tNumber of ambiguous reads\t% of ambiguous reads\tNumber of multiple alignments\t% of reads on forward strand\t% of reads on reverse strand" >> ${STATISTICS_FILE}
else
    echo -ne "Density in exons (Tags/Kb)\tDensity in introns (Tags/Kb)\tDensity in TSS up 10kb / TES down 10kb (Tags/Kb)\tNumber of assigned reads\t% of assigned reads\tNumber of no feature reads\t% of no feature reads\tNumber of ambiguous reads\t% of ambiguous reads\tNumber of multiple alignments" >> ${STATISTICS_FILE}
fi
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# launch script to do the verification of transcriptome index existence and create it if it does not exist
# => uniquely if tophat2 is used
# run with slurm if USE_SLURM=true in the config file.
#------------------------------------------------------------------------------------------------------------------

if echo ${TOPHAT_VERSION} | grep -q "tophat-2"
then
    sh $SCRIPTSDIR/src/Alignment/transcriptomeIndexVerification.sh -c ${CONFIG_FILE} -o ${OUTDIR} > ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX
fi
wait

# when the transcriptome-index verification/creation is finished, the corresponding log file is moved in the "Logs" directory.
mv ${OUTDIR}/NOHUP_TRANSCRIPTOME-INDEX ${OUTDIR}/Logs/NOHUP_TRANSCRIPTOME-INDEX



###################################################################################################################
# Launch of 2 scripts / sample
# 1) if paired = no, the RNA-seq experiment has one fastq file per sample so => one QualityAnalysis is performed
# 	 and one file is used in RNAseqPipeline1sample.
# 2) if paired = yes, the RNA-seq experiment has two fastq files per sample (.R1 & .R2) so, two quality analysis
# 	 are perfomed (one per file) and the two files are used in RNAseqPipeline1sample (but one
# 	 file is given in parameter in order to recover the common samplename, without R1 & R2).
# 3) if spikes = yes, RNAseqPipeline1sample_WithoutSpikes.sh is launch instead of RNAseqPipeline1sample & 
#	 QualityAnalysis.sh. This script separates reads of spikes from other reads and then launch 
# 	 RNAseqPipeline1sample & QualityAnalysis.sh
#------------------------------------------------------------------------------------------------------------------

########################
# START OF THE FOR LOOP
#-----------------------

if [ "$PAIRED" = "no" ]
then
	echo "***********************************************************************************"
	echo "`date`: Start analysis of all single-end samples from $DATADIR"
	# all the data files are analyzed
	for file in `ls $DATADIR/*.fastq.gz`
	do

		name=$( get_name ${file} )
		samplename=$( get_sample_name ${file} )

		echo "`date` : Processing of $samplename with the path ${name}.fastq.gz"
		
		if [ "$SPIKES" = "yes" ]
		then
		
			# when we have spikes, we process different steps to remove them from the raw fastq files
			# and for that, we need a fastq file which will contain the data without spike's reads, here:
			file_without_spikes_in_GZ_format="${OUTDIR}/Processed_data/${samplename}_without_Spikes.fastq.gz"
			#we recover the name and samplename of these files, as previously
			name_file_without_spikes=$( get_name ${file_without_spikes_in_GZ_format} )
			#sample name
			samplename_file_without_spikes=$( get_sample_name ${file_without_spikes_in_GZ_format} )

			echo "`date` : Because we have Spikes, we create file which will contains all data without spike's reads.
			So, processing of $samplename_file_without_spikes with the path ${name_file_without_spikes}.fastq.gz"

			# the script RNAseqPipeline1Sample_Without_Spikes performs the two steps performed if there are no spikes : QualityAnalysis & RNAseqPipeline1sample
			# but there is a supplementary step where reads corresponding to spikes are separated in their own file and all reads different of spikes are put in the file_without_spikes_in_GZ_format.
			cmd="sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample_WithoutSpikes.sh --raw-fastq ${file} --processed-fastq ${file_without_spikes_in_GZ_format} --raw-name ${name} --raw-samplename ${samplename} --config-file ${CONFIG_FILE} --outdir ${OUTDIR}"

			echo "#------------------------------------------"
			echo "Analysis of all samples without Spikes on $samplename_file_without_spikes :"
			echo "$cmd"

			ids[${#ids[*]}]=$(run -n "RNAseqPipeline1Sample_Without_Spikes" -p "$SLURM_PARTITION" -l "${OUTDIR}/Logs/nohup.$samplename_file_without_spikes" -w "$OUTDIR" $cmd)

		else
		
			echo "`date` : There are no Spikes in this experiment, continue with the ${name}.fastq.gz file."
			# QualityAnalysis.sh
			echo "#------------------------------------------"
			echo "Quality analysis command on $samplename :"
			echo "bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file $file --config-file ${CONFIG_FILE} --outdir ${OUTDIR} > ${OUTDIR}/Logs/nohup.$samplename &"
			echo "#------------------------------------------"

			bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file $file --config-file ${CONFIG_FILE} --outdir ${OUTDIR} > ${OUTDIR}/Logs/nohup.$samplename &

			# RNAseqPipeline1sample.sh
			# launch with slurm
			cmd="sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh --fastq-file ${file} --config-file ${CONFIG_FILE} --outdir $OUTDIR"
			
			echo "#------------------------------------------"
			echo "RNAseqPipeline1sample on $samplename :"
			echo "$cmd"
			echo "#------------------------------------------"

			ids[${#ids[*]}]=$(run -n "RNAseqPipeline1sample_${samplename}" -p "$SLURM_PARTITION" -l "${OUTDIR}/Logs/nohup.$samplename" -w "${OUTDIR}" $cmd)

		fi
	
	done

	# wait that all the tmp files are created by slurm and wait that all jobs end (when all the tmp files are removed)
	wait_until_jobs_end ${OUTDIR}
	wait

	echo " `date`: alignment of spikes and separation in two FASTQ files : DONE "

elif [ "$PAIRED" = "yes" ]
then

	echo "***********************************************************************************"
	echo "`date`: Start analysis of all paired-end samples from $DATADIR"
	# all the data files (with R1 suffix) are scanned
	# as both files correspond to one sample, if we consider all the data files
	# we analyze twice too many samples
	for file in `ls $DATADIR/*.R1.fastq.gz`
	do
	
		#path/to/sample_name (= file name without final .R1.fastq.gz)
		name=$( get_name ${file} | sed -e 's/\.R1//' )
		#sample name (the common samplename of two fastq files, without R1 or R2)
		samplename=$( get_sample_name ${file} | sed -e 's/\.R1//' )
		
		if [ "$SPIKES" = "yes" ]
		then
		
			#path/to/file without spike's reads
			file_without_spikes_in_GZ_format="${OUTDIR}/Processed_data/${samplename}_without_Spikes.fastq.gz"
			name_file_without_spikes=$( get_name ${file_without_spikes_in_GZ_format} )
			#sample name
			samplename_file_without_spikes=$( get_sample_name ${file_without_spikes_in_GZ_format} )

			cmd="sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample_WithoutSpikes.sh --raw-fastq ${file} --processed-fastq ${file_without_spikes_in_GZ_format} --raw-name ${name} --raw-samplename ${samplename} --config-file ${CONFIG_FILE} --outdir ${OUTDIR}"
			
			echo "#------------------------------------------"
			echo "Analysis of all samples without Spikes on $samplename_file_without_spikes :"
			echo "$cmd"
			echo "#------------------------------------------"

			ids[${#ids[*]}]=$(run -n "Analysis_All_Samples_Without_Spikes_${samplename}" -p "$SLURM_PARTITION" -l "${OUTDIR}/Logs/nohup.$samplename_file_without_spikes" -w "$OUTDIR" $cmd)

		else

			# RNAseqPipeline1sample.sh is performed with the common file name of paired files
			# to make the alignment with the two files and to have a common .sam file (with the two reads)
			cmd="sh $SCRIPTSDIR/src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh --fastq-file ${name}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir $OUTDIR"
			
			echo "#------------------------------------------"
			echo "RNAseqPipeline1sample on $samplename :"
			echo "$cmd"
			echo "#------------------------------------------"

			ids[${#ids[*]}]=$(run -n "RNAseqPipeline1sample_${samplename}" -p "$SLURM_PARTITION" -l "${OUTDIR}/Logs/nohup.${samplename}" -w "${OUTDIR}" $cmd)

			# 2 QualityAnalysis.sh because there are two fastq files in paired-end, so we want to know the quality of each read (R1 and R2)
			# QualityAnalysis.sh on .R1 file	
			echo "#------------------------------------------"
			echo "Quality analysis command on $samplename.R1 :"
			echo "bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} > ${OUTDIR}/Logs/nohup.$samplename"
			echo "#------------------------------------------"

			bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name}.R1.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} > ${OUTDIR}/Logs/nohup.$samplename

			wait

			# QualityAnalysis.sh on .R2 file		
			echo "#------------------------------------------"
			echo "Quality analysis command on $samplename.R2 :"
			echo "bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name}.R2.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.$samplename &"
			echo "#------------------------------------------"
			bash $SCRIPTSDIR/src/Quality/QualityAnalysis.sh --fastq-file ${name}.R2.fastq.gz --config-file ${CONFIG_FILE} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.$samplename &

		fi
		
		echo
		
	done

	# wait that all the tmp files are created by slurm and wait that all jobs end (when all the tmp files are removed)
	wait_until_jobs_end ${OUTDIR}
	wait
	
	echo " `date`: alignment of spikes and separation in two FASTQ files : DONE "
	
else

	print_err_and_exit "The paired-end variable is misnamed - EXIT"
fi

#####################
# END OF THE FOR LOOP
#--------------------
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# if spikes = yes, the .sam file creates by RNAseqPipeline1sample_WithoutSpikes.sh and containing all the reads
# of spikes is used to generates some quality analyses.
#------------------------------------------------------------------------------------------------------------------

if [ "$SPIKES" = "yes" ]
then
	
	for samfile in `ls ${OUTDIR}/Alignment/*/*_Spikes.sam`
	do

		if grep -q "@SQ" ${samfile}
		then
			echo " ${samfile} header OK "
		else
			samplename=$( basename $samfile | sed "s/_Spikes\.sam//g" )
			# creation of bam file from sam file
			SAMFILEOK=${TMP_DIR}tmp1.txt

			cat ${SCRIPTSDIR_SPIKES}/Files/SAM_header_spikes.txt ${samfile} > ${SAMFILEOK}
			wait
			mv ${SAMFILEOK} ${samfile}
		fi

		echo "$BIN/$SAMTOOLS_VERSION/samtools view -bt ${SPIKES_PATH}/Fasta/ERCC92.fa ${samfile} > ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam"
		$BIN/$SAMTOOLS_VERSION/samtools view -bt ${SPIKES_PATH}/Fasta/ERCC92.fa ${samfile} > ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam
		wait

		echo "$BIN/$SAMTOOLS_VERSION/samtools sort ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes_sorted"
		$BIN/$SAMTOOLS_VERSION/samtools sort ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes_sorted
		
		cp ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes_sorted.bam ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam
		wait
		# create de .bai file
		echo "$BIN/$SAMTOOLS_VERSION/samtools index ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam"
		$BIN/$SAMTOOLS_VERSION/samtools index ${OUTDIR}/Alignment/${samplename}_without_Spikes/${samplename}_Spikes.bam
		
	done
	wait
				
	echo "#------------------------------------------"
	echo "Spikes analysis on $samplename :"
	echo "sh $SCRIPTSDIR/src/Spikes/SpikesAllSamplesAnalysis.sh -c $CONFIG_FILE -o $OUTDIR >> ${LOGFILE_SPIKES} &"
	echo "#------------------------------------------"

	sh $SCRIPTSDIR/src/Spikes/SpikesAllSamplesAnalysis.sh -c $CONFIG_FILE -o $OUTDIR >> ${LOGFILE_SPIKES} &

	wait

fi
#------------------------------------------------------------------------------------------------------------------

echo "Analysis of all samples from $DATADIR complete."
echo "***********************************************************************************"



###################################################################################################################
# Generates a tab-delimited text file containing the statistics for all samples
# concatenate all stat files and generate corresponding tab delimited text file
# there is one summary table for the quality assessment : stats_summary.txt
# and one summary table for the reads alignment : tophat_stats_summary.txt
#------------------------------------------------------------------------------------------------------------------

echo "Creation of tab-delimited text files with the statistics on all sequences and on all tophat results."
# #statistics on the sequences
cd $OUTDIR/Quality || echo "Cannot go to $OUTDIR/Quality - Folder doesn't exist ?"
cat stats_*.txt > stats.txt
echo "source(\"$SCRIPTSDIR/src/Quality/QualityInfos.R\");dir = \"$OUTDIR/Quality\";spikes = \"$SPIKES\";qualityinfos(dir, spikes) " | $RBIN --no-save --no-restore --quiet
# #statistics on tophat results
cd $OUTDIR/Alignment || echo "Cannot go to $OUTDIR/Alignment - Folder doesn't exist ?"
cat tophat_stats_*.txt > tophat_stats.txt

if echo ${TOPHAT_VERSION} | grep -q "tophat-2"
then
	tophat2V="yes"
else

	tophat2V="no"

fi
echo "source(\"$SCRIPTSDIR/src/Alignment/TophatInfos.R\");dir = \"$OUTDIR/Alignment\";paired = \"$PAIRED\";tophat2 = \"$tophat2V\";tophatinfos(dir,paired,tophat2) " | $RBIN --no-save --no-restore --quiet

wait
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# Analysis of HTseq results and concatenation of all result files
#------------------------------------------------------------------------------------------------------------------

echo "#------------------------------------------"
echo "Analysis of HTSeq results and concatenation of all HTSeq result files :"
echo "sh $SCRIPTSDIR/src/Quantification/HTseqResultsAnalysis.sh ${OUTDIR}"
echo "#------------------------------------------"

sh $SCRIPTSDIR/src/Quantification/HTseqResultsAnalysis.sh ${OUTDIR}

wait

sort_design_file $DESIGN_FILE_REPORT $OUTDIR || { echo "/!\ /!\ sort_design_file failed - be careful that the order in htseq_combined is the same than design_file when analyzing /!\ /!\ RNAseqPipeline continue ..."; }

test_empty_file $DESIGN_FILE_REPORT
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Generation of overall automatic report 
#------------------------------------------------------------------------------------------------------------------
# create log file of automatic report
touch ${LOGFILE_REPORT}

echo "#------------------------------------------"
echo "Generation of overall automatic report :"
echo "sh ${REPORTDIR}/RNAseqReport.sh -c $CONFIG_FILE -o $OUTDIR -d $DATADIR -f $DESIGN_FILE_REPORT &>> ${LOGFILE_REPORT}"
echo "#------------------------------------------"

sh ${REPORTDIR}/RNAseqReport.sh -c $CONFIG_FILE -o $OUTDIR -d $DATADIR -f $DESIGN_FILE_REPORT &>> ${LOGFILE_REPORT}
wait

cd ${OUTDIR}/Report

#Generation of automatic report (pdflatex) and the corresponding bibliography (bibtex)
echo "#------------------------------------------"
echo "Compilation of overall latex reports :"
echo "compile_latex_report ${RBIN} ${OUTDIR} ${PROJECT_NUMBER} &>> ${LOGFILE_REPORT}"
echo "#------------------------------------------"

compile_latex_report ${RBIN} ${OUTDIR} ${PROJECT_NUMBER} &>> ${LOGFILE_REPORT}

wait

echo "`date`: generation of auto report : done"
#------------------------------------------------------------------------------------------------------------------



####################################################################################################################
# Generation of normalized wig & tdf files after recovering of sizeFactor + generation of the corresponding .bw file
#-------------------------------------------------------------------------------------------------------------------

touch ${LOGFILE_WIG_TDF}

echo "`date`: Wig & tdf files generation"

# if paired=no, all the data files are scanned to generate one wig file and one tdf file per sample
if [ "$PAIRED" = "no" ]
then

	if [ "$SPIKES" = "yes" ]
	then
	
		listFiles=`ls $OUTDIR/Processed_data/*_without_Spikes.fastq.gz`
	
	else
	
		listFiles=`ls $DATADIR/*.fastq.gz`
	
	fi
	
# if paired=yes, only the .R1 data files are scanned in order to recover the common filename and
# to generate one wig file and one tdf file per sample
elif [ "$PAIRED" = "yes" ]
then
	if [ "$SPIKES" = "yes" ]
	then
		listFiles=`ls $OUTDIR/Processed_data/*_without_Spikes.R1.fastq.gz`
	else
		listFiles=`ls $DATADIR/*.R1.fastq.gz`
	fi
	
else

	print_err_and_exit "Problem with files list for wig files generation"
	
fi


for file in ${listFiles}
do
	cmd="source $SCRIPTSDIR/src/Utilities/utils.sh && wigTdfGeneration ${file} ${CONFIG_FILE} ${OUTDIR};"
	ids[${#ids[*]}]=$(run -n "Wig_generation_${file}" -p "$SLURM_PARTITION" -l "$LOGFILE_WIG_TDF" -w "${OUTDIR}/Alignment" $cmd)
done

# wait that all the tmp files are created and wait that all jobs end
wait_until_jobs_end ${OUTDIR}
wait
#-------------------------------------------------------------------------------------------------------------------



####################################################################################################################
# Generation of gene body coverage graph
#-------------------------------------------------------------------------------------------------------------------

sh $SCRIPTSDIR/src/Quality/GeneBodyCoverage.sh --config-file ${CONFIG_FILE} --outdir ${OUTDIR}

cp $DESIGN_FILE_REPORT $DESIGN_FILE

rm -R $OUTDIR/wait_run_end/

#-------------------------------------------------------------------------------------------------------------------

echo "#------------------------------------------"
echo "`date`"
echo "RNAseqPipeline finished successfully !"
echo "#------------------------------------------"
exit 0


