#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
#This script analyze one gzipped fastq file containing RNAseq reads :
#1. mapping with Tophat and generation of a file containing statistics on tophat result
#2. quantification with HTSeq and generation of a file containing statistics on HTseq result
#3. different quality assessment in function of the parameters (stranded protocol, paired-end reads, ...)
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh RNAseqPipeline1sample.sh -f fastq-file -c config_file -o outdir [-h]

		Description:
		    This script analyze one gzipped fastq file containing RNAseq reads :
			1. mapping with Tophat and generation of a file containing statistics on tophat result
			2. quantification with HTSeq and generation of a file containing statistics on HTseq result

		Options:
		 	-f, --fastq-file a fastq file in .gz format, located in the rawdata folder (or Processed_data if Spikes=yes)
		        This option is required. 
		        Generally of type : /path/to/Num_project_Name_author/rawdata(or Processed_data)
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
ARGS=$(getopt -o "f:c:o:h" --long "fastq-file:,config-file:,outdir:,help" -- "$@" 2> /dev/null)
[ $? -ne 0 ] && \
	echo "Error in the argument list." "Use -h or --help to display the help." >&2 && \
	exit 1
eval set -- "$ARGS"
while true
do
	case "$1" in
		-f|--fastq-file)
			FILE="$2"
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



###################################################################################################################
# Checking the input parameter
#------------------------------------------------------------------------------------------------------------------
[ "$FILE" == "" ] || [ "$CONFIG_FILE" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --fastq-file, --config-file and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Source the other parameters
#------------------------------------------------------------------------------------------------------------------
source $CONFIG_FILE
source "$SCRIPTSDIR/src/Utilities/utils.sh"
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# Start of 1 sample analysis
#------------------------------------------------------------------------------------------------------------------
echo "***********************************************************************************"
echo "Analysis of file $FILE"
echo
#path/to/sample_name (= file name without final .fastq.gz)
name=$( get_name ${FILE} | sed -e 's/\.R1//' )
#sample name
samplename=$( get_sample_name ${FILE} | sed -e 's/\.R1//' )
# to wait until the script ends.
touch $OUTDIR/wait_run_end/${samplename}_tmp.txt
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# TopHat or TopHat2 alignment
# Launch TopHat on one fastq sequence file (or two fastq files if paired=yes)
# TopHat is launched with options -p nb_of_processors -G gtf_file
# results are stored in $outdir/Alignment/*/*
# where * correspond to the sample name (from *.fastq.gz)
#------------------------------------------------------------------------------------------------------------------

echo "****************************************"
echo "`date`: Tophat analysis"
echo "****************************************"

cd $OUTDIR/Alignment

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

# path to bowtie and samtools for TopHat
export PATH=$PATH:$BIN/$BOWTIE_VERSION/:$BIN/$SAMTOOLS_VERSION/

# choosing the tophat option depending on whether the protocol is stranded or not
if [ "$STRANDED" = "reverse" ]
then
    librarytype="fr-firststrand"
elif [ "$STRANDED" = "no" ]
then
    librarytype="fr-unstranded"
else
    print_err_and_exit "The STRANDED variable is misnamed. Exit."
fi

# if TopHat1 is used, no transcriptome-index are used and the alignment time is longer
if echo ${TOPHAT_VERSION} | grep -q "tophat-1"
then
	
	# if the RNA-seq libraries are not paired-end, there is only one fastq file for TopHat alignment
	if [ "$PAIRED" = "no" ]
	then
		echo "#------------------------------------------------------------------"
		echo "$TOPHAT_VERSION command : Alignment of reads on reference genome =>"
		echo "$BIN/$TOPHAT_VERSION/tophat --library-type=$librarytype -p $NBPROC --GTF=${GTF_FILE} -o $samplename ${BOWTIE_INDEXES}/${GENOME_VERSION} $FILE 2>> ${OUTDIR}/Logs/nohup.${samplename}"
		echo "#------------------------------------------------------------------"
		$BIN/$TOPHAT_VERSION/tophat --library-type=$librarytype -p $NBPROC --GTF=${GTF_FILE} -o $samplename ${BOWTIE_INDEXES}/${GENOME_VERSION} $FILE 2>> ${OUTDIR}/Logs/nohup.${samplename}
	
	# if the RNA-seq libraries are paired-end, there are two fastq files for TopHat alignment
	elif [ "$PAIRED" = "yes" ]
	then
		echo "#------------------------------------------------------------------"
		echo "$TOPHAT_VERSION command : Alignment of reads on reference genome =>"
		echo "$BIN/$TOPHAT_VERSION/tophat --library-type=$librarytype --mate-inner-dist=$FRAGMENT_SIZE -p $NBPROC --GTF=${GTF_FILE} -o $samplename ${BOWTIE_INDEXES}/${GENOME_VERSION} ${name}.R1.fastq.gz ${name}.R2.fastq.gz 2>> ${OUTDIR}/Logs/nohup.${samplename}"
		echo "#------------------------------------------------------------------"
		$BIN/$TOPHAT_VERSION/tophat --library-type=$librarytype --mate-inner-dist=$FRAGMENT_SIZE -p $NBPROC --GTF=${GTF_FILE} -o $samplename ${BOWTIE_INDEXES}/${GENOME_VERSION} ${name}.R1.fastq.gz ${name}.R2.fastq.gz 2>> ${OUTDIR}/Logs/nohup.${samplename}
	
	else

		print_err_and_exit "The paired variable is misnamed. Exit."

	fi
	
# if TopHat2 is used, transcriptome-index are used and the alignment time is faster (and the alignment better !)
elif echo ${TOPHAT_VERSION} | grep -q "tophat-2"
then
	
	if [ "$REALIGN" = "no" ]
	then

		# Launch TopHat2
		# --no-coverage-search option is used because coverage-search takes too much time or memory
		# and is necessary for short reads (< 40-bp) which is not our case.
		# --transcriptome-index option is used to reduce the time of alignment more than 30 minutes.
		# it is possible to use the --read-realign-edit-dist option, see config.sh file

		# if the RNA-seq libraries are not paired-end, there is only one fastq file for TopHat2 alignment
		if [ "$PAIRED" = "no" ]
		then

			echo "#------------------------------------------------------------------"
			echo "$TOPHAT_VERSION command : Alignment of reads on reference genome =>"
			echo "$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype -p $NBPROC --no-coverage-search --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} $FILE 2>> ${OUTDIR}/Logs/nohup.${samplename}"
			echo "#------------------------------------------------------------------"
			# possible de rajouter l'option --read-realign-edit-dist=0
			$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype -p $NBPROC --no-coverage-search --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} $FILE 2>> ${OUTDIR}/Logs/nohup.${samplename}
		
		# if the RNA-seq libraries are paired-end, there are two fastq files for TopHat2 alignment
		elif [ "$PAIRED" = "yes" ]
		then
		
			echo "#------------------------------------------------------------------"
			echo "$TOPHAT_VERSION command : Alignment of reads on reference genome =>"
			echo "$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype --mate-inner-dist=$FRAGMENT_SIZE -p $NBPROC --no-coverage-search --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} ${name}.R1.fastq.gz ${name}.R2.fastq.gz 2>> ${OUTDIR}/Logs/nohup.${samplename}"
			echo "#------------------------------------------------------------------"
			$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype --mate-inner-dist=$FRAGMENT_SIZE -p $NBPROC --no-coverage-search --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} ${name}.R1.fastq.gz ${name}.R2.fastq.gz 2>> ${OUTDIR}/Logs/nohup.${samplename}
		
		else

			print_err_and_exit "The paired variable is misnamed. Exit."

		fi

	elif [ "$REALIGN" = "yes" ]
	then

		# Launch TopHat2
		# --no-coverage-search option is used because coverage-search takes too much time or memory
		# and is necessary for short reads (< 40-bp) which is not our case.
		# --transcriptome-index option is used to reduce the time of alignment more than 30 minutes.
		# it is possible to use the --read-realign-edit-dist option, see config.sh file

		# if the RNA-seq libraries are not paired-end, there is only a fastq file for TopHat2 alignment
		if [ "$PAIRED" = "no" ]
		then
			echo "#------------------------------------------------------------------"
			echo "$TOPHAT_VERSION command : Alignment of reads on reference genome =>"
			echo "$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype -p $NBPROC --no-coverage-search --read-realign-edit-dist=0 --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} $FILE 2>> ${OUTDIR}/Logs/nohup.${samplename}"
			echo "#------------------------------------------------------------------"
			# possible de rajouter l'option --read-realign-edit-dist=0
			$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype -p $NBPROC --no-coverage-search --read-realign-edit-dist=0 --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} $FILE 2>> ${OUTDIR}/Logs/nohup.${samplename}

		# if the RNA-seq libraries are paired-end, there are two fastq files for TopHat2 alignment
		elif [ "$PAIRED" = "yes" ]
		then

			echo "#------------------------------------------------------------------"
			echo "$TOPHAT_VERSION command : Alignment of reads on reference genome =>"
			echo "$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype --mate-inner-dist=$FRAGMENT_SIZE -p $NBPROC --no-coverage-search --read-realign-edit-dist=0 --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} ${name}.R1.fastq.gz ${name}.R2.fastq.gz 2>> ${OUTDIR}/Logs/nohup.${samplename}"
			echo "#------------------------------------------------------------------"
			$BIN/$TOPHAT_VERSION/tophat2 -o $samplename --library-type=$librarytype --mate-inner-dist=$FRAGMENT_SIZE -p $NBPROC --no-coverage-search --read-realign-edit-dist=0 --transcriptome-index=${TRANSCRIPTOME_INDEX} ${BOWTIE_INDEXES}/${GENOME_VERSION} ${name}.R1.fastq.gz ${name}.R2.fastq.gz 2>> ${OUTDIR}/Logs/nohup.${samplename}
		
		else

			print_err_and_exit "The paired variable is misnamed. Exit."

		fi


	fi

else

	print_err_and_exit "Tophat version is not known. Exit."

fi

wait
wait

#Analysis of TopHat results : FINISHED
echo "Analysis of TopHat results."

#directory where TopHat results are located for this sample
dir=${OUTDIR}/Alignment/${samplename}
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# Analysis of TopHat results.
# calculation of some statistics as total number of reads, total number of unique (or multiple) aligned reads ...
#------------------------------------------------------------------------------------------------------------------
echo -n "Directory | " >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt
echo $dir >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt

####### TOTAL NUMBER OF READS #######

if [ "$PAIRED" = "yes" ]
then
	if [ "$SPIKES" = "yes" ]
	then
		# if spikes = yes, we calculate the total number of PAIRS (because paired-end sequencing) AND 
		# the total number of PAIRS without Spikes to determine the number of spikes' pairs
		raw_samplename=$( get_sample_name ${FILE} | sed -e 's/\.R1//' | sed -e 's/_without_Spikes//' )
		totalNumberReads=$( zgrep -c "@HWI" "$OUTDIR/../${DATA_FOLDER}/${raw_samplename}.R1.fastq.gz" )
		totalNumberReadsWithoutSpikes=$( zgrep -c "@HWI" "$OUTDIR/Processed_data/${samplename}.R1.fastq.gz" )
		filteredOut=$( head -4 $OUTDIR/Alignment/${samplename}/logs/prep_reads.log | tail -n 2 | cut -d" " -f1 | awk '{ SUM += $1} END { printf("%d\n",SUM/NR) }' )

		echo -n "Nb of reads | ${filteredOut} out of ${totalNumberReadsWithoutSpikes} reads have been filtered out
" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt
	else
		# else we calculate only the total number of PAIRS (pairs because we calculate only on .R1 file)
		totalNumberReads=$( zgrep -c "@HWI" "$OUTDIR/../${DATA_FOLDER}/${samplename}.R1.fastq.gz" )
		filteredOut=$( head -4 $OUTDIR/Alignment/${samplename}/logs/prep_reads.log | tail -n 1 | cut -d" " -f1 )
		#nb of reads filtered out by tophat
		echo -n "Nb of reads | ${filteredOut} out of ${totalNumberReads} reads have been filtered out
" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt
	fi
else
	if [ "$SPIKES" = "yes" ]
	then
		# if spikes = yes, we calculate the total number of reads AND the total number of reads without Spikes to determine the number of spikes' reads
		raw_samplename=$( get_sample_name ${FILE} | sed -e 's/\.R1//' | sed -e 's/_without_Spikes//' )
		totalNumberReads=$( zgrep -c "@HWI" "$OUTDIR/../${DATA_FOLDER}/${raw_samplename}.fastq.gz" )
		totalNumberReadsWithoutSpikes=$( zgrep -c "@HWI" "$OUTDIR/Processed_data/${samplename}.fastq.gz" )
		filteredOut=$( head -3 $OUTDIR/Alignment/${samplename}/logs/prep_reads.log | tail -n 1| cut -d" " -f1 )

		echo -n "Nb of reads | ${filteredOut} out of ${totalNumberReadsWithoutSpikes} reads have been filtered out
" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt
	else
		# else we calculate only the total number of reads
		totalNumberReads=$( zgrep -c "@HWI" "$OUTDIR/../${DATA_FOLDER}/${samplename}.fastq.gz" )
		filteredOut=$( head -3 $OUTDIR/Alignment/${samplename}/logs/prep_reads.log | tail -n 1| cut -d" " -f1 )
		#nb of reads filtered out by tophat
		echo -n "Nb of reads | ${filteredOut} out of ${totalNumberReads} reads have been filtered out
" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt
	fi
fi

#conversion of bam to sam => to generate the alignment statistics
$BIN/$SAMTOOLS_VERSION/samtools view $dir/accepted_hits.bam -o $dir/accepted_hits.sam 
wait

####### STATISTICS ON ALIGNMENT #######
# calculation of statistics on alignment file (.sam) : eg nb of aligned reads ...
# with TopHat1, we must do the calculations ourselves.
# with TopHat2, the "align_summary.txt" file is generated for each sample and contain the different statistics, so we take the statistics in this file.
# but with TopHat2, the statistics for paired-end in the "align_summary.txt" file is not well documented, so for paired-end we calculate also the statistics ourselves.

if echo ${TOPHAT_VERSION} | grep -q "tophat-1"
then
			
	#count nb of aligned reads
	echo -n "Nb of aligned reads | " >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt #(note that 2 reads could have the same sequence)
	cut -f 1 $dir/accepted_hits.sam | sort -T $TMP_DIR | uniq | wc -l >> $OUTDIR/Alignment/tophat_stats_$samplename.txt
		
	#count nb of uniquely aligned reads
	echo -n "Nb of uniquely aligned reads | " >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt #(note that 2 reads could have the same sequence)
	cut -f 1 $dir/accepted_hits.sam | sort -T $TMP_DIR | uniq -u | wc -l >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt 
		
	#count nb of multiple aligned reads
	#caution, all reads with more than 40 alignments are not reported by tophat
	echo -n "Nb of multiple aligned reads | " >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt #(note that 2 reads could have the same sequence)
	cut -f 1 $dir/accepted_hits.sam | sort -T $TMP_DIR | uniq -d | wc -l >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt 

elif echo ${TOPHAT_VERSION} | grep -q "tophat-2"
then

	if [ "$PAIRED" = "no" ]
	then

		# from the "align_summary.txt" generated by TopHat2
		Mapped=`grep "Mapped" ${dir}/align_summary.txt | sed "s/ *.Mapped   ://" | sed "s/(.*)//" | sed "s/ //g"`
		Multiple=`grep "of these:" ${dir}/align_summary.txt | sed "s/ *.of these://" | sed "s/(.*)//" | sed "s/ //g"`
		Unique=$( bc <<<"scale=3;(${Mapped}-${Multiple})" ) 
				
		#count nb of aligned reads
		echo "Nb of aligned reads (or aligned pairs for paired-end) | ${Mapped}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt #(note that 2 reads could have the same sequence)
				
		#count nb of uniquely aligned reads
		echo "Nb of uniquely aligned reads | ${Unique}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt #(note that 2 reads could have the same sequence)
				
		#count nb of multiple aligned reads
		#caution, all reads with more than 40 alignments are not reported by tophat
		echo "Nb of multiple aligned reads | ${Multiple}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt #(note that 2 reads could have the same sequence)

	elif [ "$PAIRED" = "yes" ]
	then

		# from the "align_summary.txt" generated by TopHat2
		Mapped=`grep "Aligned pairs" ${dir}/align_summary.txt | sed "s/Aligned pairs://" | sed "s/ //g"`
		Multiple=`grep "have multiple alignments$" ${dir}/align_summary.txt | sed "s/ *.of these://" | sed "s/(.*).*//" | sed "s/ //g"`
		Unique=$( bc <<<"scale=3;(${Mapped}-${Multiple})" ) 
		# supplementary statistics give by TopHat2 for the PE experiment
		Discordant=`grep "discordant" ${dir}/align_summary.txt | sed "s/ (.*//" | sed 's/ //g'`

		#count nb of aligned reads
		echo "Nb of aligned pairs calculated by TopHat2 | ${Mapped}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt
				
		#count nb of uniquely aligned reads
		echo "Nb of uniquely aligned pairs calculated by TopHat2 | ${Unique}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt
				
		#count nb of multiple aligned reads
		#caution, all reads with more than 40 alignments are not reported by tophat
		echo "Nb of multiple aligned pairs calculated by TopHat2 | ${Multiple}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt


		# statistics calculate ourselves, and not by TopHat2, because we see a difference between the statistics given by TopHat2 in align_summary.txt
		# and the statistics calculated by us

		#count nb of aligned reads
		Mapped2=$( cut -f 1 $dir/accepted_hits.sam | sort -T $TMP_DIR | uniq | wc -l )
		echo "Nb of aligned reads | $Mapped2" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt #(note that 2 reads could have the same sequence)
		
		#count nb of uniquely aligned reads
		grep "	NH:i:1$" $dir/accepted_hits.sam > $dir/tmp.sam
		grep "NH:i:1	" $dir/accepted_hits.sam >> $dir/tmp.sam
		Unique2=$( cut -f1 $dir/tmp.sam | sort -T $TMP_DIR | uniq | wc -l )
		echo "Nb of uniquely aligned reads | $Unique2" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt
		wait
		rm $dir/tmp.sam

		#count nb of multiple aligned reads
		#caution, all reads with more than 40 alignments are not reported by tophat
		multiple2=$( bc <<<"scale=0;(${Mapped2}-${Unique2})" )
		echo "Nb of multiple aligned reads | $multiple2" >> $OUTDIR/Alignment/tophat_stats_${samplename}.txt 

		#count nb of discordant paired reads
		#caution, all reads with more than 40 alignments are not reported by tophat
		echo "Nb of discordant pair alignments (only for paired-end sequencing) | ${Discordant}" >> $OUTDIR/Alignment/tophat_stats_$samplename.txt #(note that 2 reads could have the same sequence)

	else 

		print_err_and_exit "The paired variable is misnamed. Exit."

	fi

else

	echo "Tophat version is not known"

fi

#change accepted_hits.bam to X_alignment.bam where X is the name of the directory ($dir)
mv $dir/accepted_hits.bam $dir/"$samplename"_alignment.bam 

#change junctions.bed to X_junctions.bed where X is the name of the directory ($dir)
mv $dir/junctions.bed $dir/"$samplename"_junctions.bed

#change accepted_hits.sam to X_alignment.sam where X is the name of the directory ($dir)
mv $dir/accepted_hits.sam $dir/"$samplename"_alignment.sam
	
#change insertions.bed to X_insertions.bed where X is the name of the directory ($dir)
mv $dir/insertions.bed $dir/"$samplename"_insertions.bed
	
#change deletions.bed to X_deletions.bed where X is the name of the directory ($dir)
mv $dir/deletions.bed $dir/"$samplename"_deletions.bed
	
#indexation of bam file
$BIN/$SAMTOOLS_VERSION/samtools index $dir/"$samplename"_alignment.bam
wait

echo "`date`: Tophat analysis : done "
#------------------------------------------------------------------------------------------------------------------



###################################################################################################################
# HTSeq analysis
# Launch htseq-count on the sam file obtained previously by TopHat
# Results are stored in $OUTDIR/Alignment/$samplename directory
# htseq-count is launched with intersection-nonempty mode
# name of the results file : htseq_name.txt
#------------------------------------------------------------------------------------------------------------------
echo "****************************************"
echo "`date`: HTSeq analysis"
echo "****************************************"

# export of the python path where find the eggs.
export PYTHONPATH=${PYTHONPATH}/${HTSEQ_PYTHONPATH}/

samfile="$dir/${samplename}_alignment.sam"

# paired reads must be sorted by name for that HTSeq runs (ie the paired reads must be one after the other)
# this sorting is no longer necessary with the new version of htseq-count, which is ok with positions sorting
# BUT with sorting by position, all the alignment is put in a buffer ... not optimal. So, we keep the option "sort by name"
if [ "$PAIRED" = "yes" ]
then

	# sort by name with samtools sort
	$BIN/$SAMTOOLS_VERSION/samtools sort -n $dir/"$samplename"_alignment.bam $dir/"$samplename"_alignment_sort_by_name
	# create the corresponding sam
	$BIN/$SAMTOOLS_VERSION/samtools view -h -o $dir/"$samplename"_alignment_sort_by_name.sam $dir/"$samplename"_alignment_sort_by_name.bam
	wait
	rm $dir/"$samplename"_alignment_sort_by_name.bam
	samfile="$dir/${samplename}_alignment_sort_by_name.sam"

fi
wait

echo "#------------------------------------------------------------------"
echo "$HTSEQ_VERSION command : Genes' quantification =>"
echo "$BIN/$HTSEQ_VERSION/htseq-count -m intersection-nonempty -s $STRANDED ${samfile} $GTF_FILE > $dir/${samplename}_htseq.txt 2>> ${OUTDIR}/Logs/nohup.${samplename}"
echo "#------------------------------------------------------------------"

$BIN/$HTSEQ_VERSION/htseq-count -m intersection-nonempty -s $STRANDED ${samfile} $GTF_FILE > $dir/${samplename}_htseq.txt 2>> ${OUTDIR}/Logs/nohup.${samplename}

wait
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# READS DISTRIBUTION
# Analysis of reads distribution in the exonic , intronic and intergenic regions
# The script "read_distribution.py" generates a file with
# the density of tags (Tags/Kb) in each region (intron, intergenic exonic)
# Requires the existence of the gtf reference file in bed format
# If this bed file does not exist, it is automatically created with the function "create_bed_from_gtf" in utils.sh
#------------------------------------------------------------------------------------------------------------------

resultsFileRSeQC=$OUTDIR/Quality/RSeQC/resultsTable_${samplename}.txt
touch ${resultsFileRSeQC}

echo "#------------------------------------------------------------------"
echo "$RSEQC_VERSION command : Reads distribution =>"
echo "sh $SCRIPTSDIR/src/Quality/ReadsDistribution.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} &"
echo "#------------------------------------------------------------------"
sh $SCRIPTSDIR/src/Quality/ReadsDistribution.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} &
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# STRAND SPECIFICITY => ONLY FOR STRANDED PROTOCOLS
# Analysis of strand specificity in the experiment
# The script "strand_specificity.py" generates a file with
# the fraction of reads on each strand (eh the number of reads explained by "++,--" or "+-,-+" for single-read)
# Requires the existence of the gtf reference file in bed format
# If this bed file does not exist, it is automatically created with the function "create_bed_from_gtf" in utils.sh
# Requires fasta file of the genome use. For this, we use $assembly value to convert in corresponding fasta file
# the conversion file is in RNAseqReport/Files folder, and is called "genome_conversion.txt"
#------------------------------------------------------------------------------------------------------------------
if [ "$STRANDED" = "reverse" ]
then
	
	echo "#------------------------------------------------------------------"
	echo "$RSEQC_VERSION command : Strand specificity =>"
	echo "sh $SCRIPTSDIR/src/Quality/StrandSpecificity.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} &"
	echo "#------------------------------------------------------------------"

	sh $SCRIPTSDIR/src/Quality/StrandSpecificity.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} &

fi

wait
wait
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# INNER DISTANCE => ONLY FOR PAIRED-END PROTOCOLS
# This module is used to calculate the inner distance (or insert size) between two paired RNA reads. 
# The distance is the mRNA length between two paired fragments. We first determine the genomic (DNA) size
# between two paired reads: D_size = read2_start - read1_end, then
#    if two paired reads map to the same exon: inner distance = D_size
#    if two paired reads map to different exons:inner distance = D_size - intron_size
#    if two paired reads map non-exonic region (such as intron and intergenic region): inner distance = D_size
#    The inner_distance might be a negative value if two fragments were overlapped.
#------------------------------------------------------------------------------------------------------------------
if [ "$PAIRED" = "yes" ]
then

	echo "#------------------------------------------------------------------"
	echo "$RSEQC_VERSION command : Inner Distance for Paired-end data =>"
	echo "sh $SCRIPTSDIR/src/Quality/InnerDistance.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} &"
	echo "#------------------------------------------------------------------"
	
	sh $SCRIPTSDIR/src/Quality/InnerDistance.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} &

fi
#------------------------------------------------------------------------------------------------------------------

cd $OUTDIR

###################################################################################################################
# OVERALL STATISTICS FOR EACH SAMPLE
# create a table of all the statistics for each sample
# e.g sequencing quality statistics, tophat aligment statistics, reads distribution ...
# available in ${OUTDIR}/Quality/all_statistics_results.txt (statisticsFile)
#------------------------------------------------------------------------------------------------------------------

echo "#------------------------------------------------------------------"
echo "$RSEQC_VERSION command : Overall statistics for $samplename =>"
echo "sh $SCRIPTSDIR/src/Quality/OverallStatistics.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --fastq-file ${FILE} --resultsRSeQC ${resultsFileRSeQC} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} & "
echo "#------------------------------------------------------------------"

sh $SCRIPTSDIR/src/Quality/OverallStatistics.sh --config-file ${CONFIG_FILE} --samplename ${samplename} --fastq-file ${FILE} --resultsRSeQC ${resultsFileRSeQC} --outdir ${OUTDIR} >> ${OUTDIR}/Logs/nohup.${samplename} & 

wait
#------------------------------------------------------------------------------------------------------------------

# to wait until the script ends.
rm $OUTDIR/wait_run_end/${samplename}_tmp.txt

echo "#------------------------------------------------------------------"
echo "#------------------------------------------------------------------"
echo "#------------------------------------------------------------------"
echo "`date` : RNAseqPipeline1sample for ${samplename} is finished"
echo "#------------------------------------------------------------------"
echo "#------------------------------------------------------------------"
echo "#------------------------------------------------------------------"

