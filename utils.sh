#! /bin/bash
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine Velt
# Mail:    velt@igbmc.fr
# Date:    March 2015
#------------------------------------------------------------------------------------------------------------------
# This script contains some functions that will be usefull for the RNAseq pipeline.
###################################################################################################################


#------------------------------------------------------------------------------------------------------------------
# This function is used to create a list of directories.
#------------------------------------------------------------------------------------------------------------------
# INPUT:  
#   A list of directory names.
# OUTPUT:
#   Nothing.
# RETURN:
#   0 in success, otherwise exit 1.
#------------------------------------------------------------------------------------------------------------------
function make_directories
{
	for dir in $*
	do
		if [ ! -d "$dir" ]
		then
			mkdir -p "$dir" >> /dev/null 2>&1 
		fi
		if [ $? -ne 0 ]
		then
			print_err_and_exit "Cannot create the directory $dir" 
		fi
	done
	return 0
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function prints a list of strings on the stderr stream.
#------------------------------------------------------------------------------------------------------------------
# INPUT:  
#   A list of string.
# OUTPUT: 
#   Nothing.
# RETURN:
#   Exit status of the last command.
#------------------------------------------------------------------------------------------------------------------
function echoerr
{
	echo "ERROR: " $* >&2
	return $?
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function prints the error message and exits the program.
#------------------------------------------------------------------------------------------------------------------
# INPUT:  
#   A list of strings.
# OUTPUT: 
#   Nothing.
# RETURN:
#   Exit 1.
#------------------------------------------------------------------------------------------------------------------
function print_err_and_exit
{
	echoerr $*
	exit 1
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function print the current time
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   Nothing
# OUTPUT:
#   The actual time given like "janvier 23 2015 - 18:16:52"
# RETURN:
#   Exit state of the last command.
#------------------------------------------------------------------------------------------------------------------
function print_current_time
{
	date \+\%B" "\%d" "\%Y" - "\%X
	return $?
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function securizes the writing access of the file given in first parameter. 
# So only one job can write in it.
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - First parameters the file name where to write
#   - All the remaining parametrs are the text to write in the file given above
# OUTPUT:
#   Nothing
# RETURN:
#   Exit state of the last command 
#------------------------------------------------------------------------------------------------------------------
function echo_one_at_a_time
{
	file="$1"
	shift
	( flock -e 200 ; echo "$@" 1>&200 2>&1 ) 200>>"$file"
	return $?
}
#------------------------------------------------------------------------------------------------------------------
function cat_one_at_a_time
{
	file="$1"
	shift
	( flock -e 200 ; cat "$@" 1>&200 2>&1 ) 200>>"$file"
	return $?
}
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# Append command and the stdout and stderr in the log file given in first argument
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - The log file
#   - The command
# OUTPUT:
#   Nothing
# RETURN:
#   Exit state of the command.
#------------------------------------------------------------------------------------------------------------------
function append_cmd_to_log
{
	log=$1
	shift
	cmd="$@"
	tmp_log=$(mktemp --tmpdir=$(dirname "$log") --suffix=$(basename "$log"))

	echo_one_at_a_time "$tmp_log"
	echo_one_at_a_time "$tmp_log" "$cmd"
	echo_one_at_a_time "$tmp_log" "# Started at: $(print_current_time)"
	echo_one_at_a_time "$tmp_log" "#<--- Output: ---------------------------------------------------"
	
	bash -c "$cmd" >> "$tmp_log" 2>&1
	state=$?
	
	echo_one_at_a_time "$tmp_log"
	echo_one_at_a_time "$tmp_log" "#--------------------------------------------------------------->"
	echo_one_at_a_time "$tmp_log" "# Ended at: $(print_current_time)"
	echo_one_at_a_time "$tmp_log"
	
	cat_one_at_a_time "$log" "$tmp_log"
	rm "$tmp_log"
	
	return $state
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function takes several arguments as input:
#  -n) To specify the job name
#  -p) To specify the partition to be used
#  -x) A list node to be exclude from the run, separated by a comma
#  -O) The path to the file where to write the stdout 
#  -E) The path to the file where to write the stderr
#  -d) The dependency list jobs, separated by a comma
#  -w) The path to the working directory that you want to use
#  -t) The number of tasks per node you want to use
#  -m) The minimal memory you want to allocate per node
#  -l) Append the stdout and stderr on the specified file; this option
#      will disable the -O and -E options
#  Then the command that you want to launch
#------------------------------------------------------------------------------------------------------------------
# It returns the exit status of your command, and if you have used slurm it prints on the
# stdout the task id of your submitted job.
#------------------------------------------------------------------------------------------------------------------
function run
{
	# Initialize some variables
	job_name=""
	partition=""
	exclude=""
	stdout=""
	stderr=""
	dependency=""
	task_per_nodes=""
	working_directory=""
	memory=""
	cmd=""
	log=""

	# Read your options and stop when an option isn't recognized
	OPTIND=1
	while getopts "n:p:x:O:E:d:w:t:m:l:" option "$@" > /dev/null 2>&1
	do
		case "$option" in
			n)
				[ "$SBATCH_JOB_NAME" == "" ] || job_name="$SBATCH_JOB_NAME \"$OPTARG\"" ;;
			p)
				[ "$SBATCH_PARTITION" == "" ] || partition="$SBATCH_PARTITION $OPTARG" ;;
			x)
				[ "$SBATCH_EXCLUDE" == "" ] || exclude="$SBATCH_EXCLUDE $OPTARG" ;;
			O)
				[ "$SBATCH_STDOUT" == "" ] || [ "$log" != "" ] || stdout="$SBATCH_STDOUT \"$OPTARG\"" ;;
			E)
				[ "$SBATCH_STDERR" == "" ] || [ "$log" != "" ] || stderr="$SBATCH_STDERR \"$OPTARG\"" ;;
			d)
				[ "$SBATCH_DEPENDENCY" == "" ] || \
					dependency="$SBATCH_DEPENDENCY $SBATCH_DEPENDENCY_OPTION$(echo $OPTARG | sed "s/,/$SBATCH_DEPENDENCY_SEPARATOR/g")" 
				;;
			w)
				[ "$SBATCH_WORKING_DIRCTORY" == "" ] || working_directory="$SBATCH_WORKING_DIRCTORY \"$OPTARG\"" ;;
			t)
				[ "$SBATCH_TASKS_PER_NODE" == "" ] || task_per_nodes="$SBATCH_TASKS_PER_NODE $OPTARG" ;;
			m)
				[ "$SBATCH_MEMORY" == "" ] || memory="$SBATCH_MEMORY $OPTARG" ;;
			l)
				log="$OPTARG"
				[ "$SBATCH_STDOUT" == "" ] || stdout="$SBATCH_STDOUT /dev/null"
				[ "$SBATCH_STDERR" == "" ] || stderr="$SBATCH_STDERR /dev/null"
				;;
		esac
	done
	shift $((OPTIND-1))
	cmd="$@"

	# If you using a log file
	if [ "$log" != "" ]
	then
		cmd="source \"$SCRIPTSDIR/src/Utilities/utils.sh\" ; append_cmd_to_log \"$log\" '"$cmd"' ;"
	fi

	# Launch your command with slurm or not.
	if [ $USE_SLURM ]
	then
		id=$($SBATCH $job_name $partition $exclude $stdout $stderr $dependency $task_per_nodes $working_directory $memory <<-EOF
			#!/bin/bash
			bash -c "$cmd"
			EOF
			)
		out=$?
		id=$(echo $id | cut -d " " -f 4)
		echo $id
	else
		bash -c "$cmd"
		out=$?
	fi

	# Exit
	return $out
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function takes as input a directory and wait until all files are removed (= jobs ended).
#------------------------------------------------------------------------------------------------------------------

function wait_until_jobs_end {

	OUTDIR=$1
	sleep 30
	while [ `ls $OUTDIR/wait_run_end/ | wc -l` != 0 ]
	do
	  sleep 60
	done

}

#------------------------------------------------------------------------------------------------------------------
# This function takes as input a list of id and wait until job is not finished.
#------------------------------------------------------------------------------------------------------------------
# function wait_until_ids_end
# {
# 	if [ $USE_SLURM ] && [ $# -ne 0 ]
# 	then
# 		pattern=""
# 		for id in "$@"
# 		do
# 			pattern=" -e $id "$pattern
# 		done

# 		while : ; do
# 			squeue_output=$(squeue)
# 			out=$?
# 			grep_output=$( echo "${squeue_output}" | grep "JOBID" )
			
# 			if [ $out -eq 0 ] && [ ! -z "$grep_output" ] && [ $(echo $squeue_output | grep $pattern | wc -l) -eq 0 ]
# 			then
# 				break
# 			fi
			
# 			sleep 30
# 		done

# 	fi
# 	return 0
# }
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# This function prints the sample name once the prefix (path) and suffix (fastq.gz) have been removed from the given path.
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - A complete path to a fastq.gz file
# OUTPUT: 
#   The sample name without prefix and suffix
#------------------------------------------------------------------------------------------------------------------
function get_name
{
		#path/to/sample_name (= file name without final .fastq.gz)
		name=`echo ${1} | sed "s/\.fastq\.gz//"`
		echo "$name"
}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function prints the sample name once the suffix (fastq.gz) have been removed from the given path.
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   -  A complete path to a fastq.gz file
# OUTPUT: 
#   The name without suffix
#------------------------------------------------------------------------------------------------------------------
function get_sample_name
{

		#sample name
		samplename=`basename ${1} | sed "s/\.fastq\.gz//"`
		echo "$samplename"
}
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# This function prints the bed file "path/name" (and create this file if it doesn't exist)
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   -  Path to the gtf file
#	- Path to the scripts directory (where find gtf2bed.pl )
# OUTPUT: 
#   - The bedfile path
#------------------------------------------------------------------------------------------------------------------
function create_bed_from_gtf
{

	gtf_file=${1}
	scriptsdir=${2}

	bedfile=$( echo ${gtf_file}.bed | sed "s/\.gtf//" )
		
	if [ ! -e "$bedfile" ]
	then
		
		perl ${scriptsdir}/Utilities/gtf2bed.pl ${GTF_FILE} > ${bedfile}
		
	fi

	echo "${bedfile}"
}
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# This function prints the design_file so that samplenames are in the same order as in htseq_combined
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - Path to the raw design file
#	- Path to outdir
# OUTPUT: 
#   - the correct design file (which replace the raw design file)
#------------------------------------------------------------------------------------------------------------------
function sort_design_file {

	design_file=${1}
	outdir=${2}

	head -1 $design_file > $outdir/Scripts/tmp.txt

	for sample in `head -1 $outdir/Quantification/htseq_combined.txt | sed 's/GeneID//;s/Ensembl gene id//' | tr '\t' '\n' | sed "1d"`
	do

		grep "$sample" $design_file >> $outdir/Scripts/tmp.txt

	done

	mv $outdir/Scripts/tmp.txt $design_file

}
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
# This function verify that the file is not empty. If the file is empty, the pipeline exit => use for important files
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - Path to the file
#------------------------------------------------------------------------------------------------------------------
function test_empty_file {

	file=$1
	if [ -s $file ]; then
	  echo "file is non empty - OK - continue"
	else
	  echo "file is empty - pipeline STOP"
	  exit 1
	fi

}
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
# This function allows to create the wig & tdf files (normalized if the statistical analysis is choosed, 
# not normalized alternatively) for one sample together with the bigwig file
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - Path to the fastq file
#	- Path to the config file
#	- Path to outdir
# OUTPUT: 
#   - One wig file (normalized if STATISTICAL_ANALYSIS="yes" in config_file)
#	- One tdf file (normalized if STATISTICAL_ANALYSIS="yes" in config_file)
#	- One bigwig file
#------------------------------------------------------------------------------------------------------------------
function wigTdfGeneration {
	
	file=${1}
	config_file=${2}
	outdir=${3}
	source ${config_file}
	
	DESIGN_FILE="$outdir/Scripts/design_file.txt"
	
	#sample_name (= file name without final .fastq.gz)
	samplename=$( get_sample_name ${file} | sed -e 's/.R1//' )
	touch $outdir/wait_run_end/${samplename}_tmp_wig.txt
		
	# sizeFactor (calculated in RNAseqDataExploration script)
	# in the file $outdir/Report/file/numProject_sizeFactor.txt,
	# the sample name does not contain "_without_Spikes" => we must eliminate this word to grep the sizefactor of sample
		
	if [ -f ${outdir}/Report/files/${PROJECT_NUMBER}_sizeFactor.txt ]
	then
		
		sizeFactorUse=`cat ${outdir}/Report/files/${PROJECT_NUMBER}_sizeFactor.txt | grep ${samplename} | cut -f2`
		
	fi
	# it is tested whether normfactor is empty, that is to say if the file sizefactor 
	# could not be generated during the statistical analysis. in this case it puts the normfactor to 1 and 
	# generates unnormalized wig files 
		
	if [ -z ${sizeFactorUse} ]
	then
		sizeFactorUse=1
	fi
		
	real_Samplename=$( grep ${samplename} ${DESIGN_FILE} | cut -f2 )
	echo "****************************************"
	echo "Command that will be launched :"
	echo 'nohup sh ${SCRIPTSDIR}/src/Utilities/generateWig.sh ${outdir}/Alignment/${samplename}/${samplename}_alignment.sam ${sizeFactorUse} "track type=bedGraph name=\"${samplename}\" description=\"${real_Samplename}\" visibility=full alwaysZero=on windowingFunction=maximum" ${config_file}&'
	echo "****************************************"
	# wig file and tdf file generation
	sh ${SCRIPTSDIR}/src/Utilities/generateWig.sh ${outdir}/Alignment/${samplename}/${samplename}_alignment.sam ${sizeFactorUse} "track type=wiggle_0 name=\"${samplename}\" description=\"${real_Samplename}\" visibility=full alwaysZero=on windowingFunction=maximum" ${config_file} ${outdir}
	wait
		
		
	#split of wig files to upload on UCSC : one file per chromosome
	mkdir ${outdir}/Alignment/${samplename}/split_wigs
	cd ${outdir}/Alignment/${samplename}/split_wigs
	wig=`ls ${outdir}/Alignment/${samplename}/${samplename}*.wig`
	filewig=`echo ${wig} | awk -F/ '{print $NF}' | sed -e 's/.wig//'`
	awk 'BEGIN{file="'${outdir}/Alignment/${samplename}/split_wigs/${filewig}'""_to_supress"}/^variableStep/{file="'${outdir}/Alignment/${samplename}/split_wigs/${filewig}'""_"(++i)".wig"}{print > file}' ${wig}
	rm *to_supress
		
		
	# add of trackline specific to each split file
	for i in ` ls *.wig `
	do
		nameWig=`echo ${i} | cut -d"." -f1`
		sed -i "1itrack type=wiggle_0 name='${nameWig}' description='${real_Samplename}' visibility=full alwaysZero=on windowingFunction=maximum" ${i}
	done
	wait
		
	#creation of corresponding bigwig file
	export PATH=$PATH:${UTILITIESDIR}/UCSC_tools/
	cd ${UTILITIESDIR}/UCSC_tools/
	BigWig=`echo ${wig} | sed -e 's/.wig/.bw/'`
	wigToBigWig ${wig} ${SCRIPTSDIR}/src/Utilities/chr_length_${GENOME_VERSION}.txt ${BigWig}
		
	#Gzip SAM file, wig file and *.bed files
	gzip ${wig}
	gzip ${outdir}/Alignment/${samplename}/${samplename}_alignment.sam
	gzip ${outdir}/Alignment/${samplename}/*.bed
	
	rm $outdir/wait_run_end/${samplename}_tmp_wig.txt
		
}
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
# This function compiles all the latex reports of the RNAseq pipeline
#------------------------------------------------------------------------------------------------------------------
# INPUT:
#   - Path to the R binary
#	- Path to outdir
#	- Project number
# OUTPUT :
#	- one data exploration report and one statistical analysis report
#------------------------------------------------------------------------------------------------------------------
function compile_latex_report {

	RBIN=$1
	OUTDIR=$2
	PROJECT_NUMBER=$3

	cd ${OUTDIR}/Report

	${RBIN} CMD pdflatex ${OUTDIR}/Report/${PROJECT_NUMBER}_report.tex
	${RBIN} CMD bibtex ${OUTDIR}/Report/${PROJECT_NUMBER}_report
	${RBIN} CMD pdflatex ${OUTDIR}/Report/${PROJECT_NUMBER}_report.tex
	${RBIN} CMD pdflatex ${OUTDIR}/Report/${PROJECT_NUMBER}_report.tex

	${RBIN} CMD pdflatex ${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataExploration.tex
	${RBIN} CMD bibtex ${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataExploration
	${RBIN} CMD pdflatex ${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataExploration.tex
	${RBIN} CMD pdflatex ${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataExploration.tex

	mv *.tex *.aux *.log *.rnw *.bbl *.blg ${OUTDIR}/Report/files/tex_files
	rm ${OUTDIR}/samplename*.txt ${OUTDIR}/Report/*.texe

}
#------------------------------------------------------------------------------------------------------------------
