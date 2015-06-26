#! /bin/bash
#--------------------------------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#--------------------------------------------------------------------------------------------------------------------------------------------
# This script analyze one gzipped fastq file containing RNAseq reads (without spikes) :
# it performs quality analysis with FastQC, generation of a png file with quality graphs and a file containing basic statistics on the sample
#--------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh QualityAnalysis.sh -f fastq-file -c config_file -o outdir [-h]

		Description:
		    This script analyze one gzipped fastq file containing RNAseq reads (without spikes) :
		    it performs quality analysis with FastQC, generation of a png file with quality graphs and a file 
		    containing basic statistics on the sample

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

#------------------------------------------------------------------------------------------------------------------
# Checking the input parameter
#------------------------------------------------------------------------------------------------------------------
[ "$FILE" == "" ] || [ "$CONFIG_FILE" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --fastq-file, --config-file and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"


#------------------------------------------------------------------------------------------------------------------
#Quality analysis
#------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#- Calculate total and unique number of sequences (in $outdir/Quality/stats_$samplename.txt file)
#- Launch fastqc (in $outdir/Quality/fastqc)
#- Generate a png file $outdir/Quality with 2 graphs (in $outdir/Quality/) :
#1. boxplot of quality score in each sequencing cycle
#2. % of each nt in each sequencing cycle
#-------------------------------------------------------------------------------------------------
# for verification
echo "file : $FILE"

cd $OUTDIR/Quality

# unique sample_name (with R1 or R2 if paired-end data) (= file name without final .fastq.gz)
samplename=$( get_sample_name ${FILE} )

echo "#----------------------------------------"
echo "`date`: Quality analysis"
echo "#----------------------------------------"

echo "File |" $FILE > $OUTDIR/Quality/stats_${samplename}.txt

#name of the sequencer, eg ILLUMINA-199BA2
seqname=`zcat $FILE | head -1 | awk -F ":" '{print $1}' | awk -F "@" '{print $2}'`
echo "Sequencer |" $seqname >> $OUTDIR/Quality/stats_${samplename}.txt

#flow cell number
flowcell=`zcat $FILE | head -1 | awk -F ":" '{print $3}'`
echo "Flow cell |" $flowcell >> $OUTDIR/Quality/stats_${samplename}.txt

#line number
line=`zcat $FILE | head -1 | awk -F ":" '{print $4}'`
echo "Line |" $line >> $OUTDIR/Quality/stats_${samplename}.txt

#read length
seqlengthplus1=`zcat $FILE | head -2 | tail -1 | wc | awk -F " " '{print $3}'`
seqlength=$(( ${seqlengthplus1} - 1 ))
echo "Read length | " $seqlength >> $OUTDIR/Quality/stats_${samplename}.txt

if [ "$SPIKES" = "yes" ]
then

	#calculate total number of sequences with spikes
	raw_samplename=$( echo ${samplename} | sed "s/_without_Spikes//" )
	seqname=`zcat $OUTDIR/../${DATA_FOLDER}/${raw_samplename}.fastq.gz | head -1 | awk -F ":" '{print $1}' | awk -F "@" '{print $2}'`
	echo -n "Total number of sequences | " >> $OUTDIR/Quality/stats_${samplename}.txt
	zgrep -c "^@$seqname" $OUTDIR/../${DATA_FOLDER}/${raw_samplename}.fastq.gz >> $OUTDIR/Quality/stats_${samplename}.txt

	#calculate total number of sequences without spikes
	echo -n "Number of sequences without spikes | " >> $OUTDIR/Quality/stats_${samplename}.txt
	zgrep -c "^@$seqname" $FILE >> $OUTDIR/Quality/stats_${samplename}.txt

else

	#calculate total number of sequences without spikes
	echo -n "Total number of sequences | " >> stats_${samplename}.txt
	zgrep -c "^@$seqname" $FILE >> $OUTDIR/Quality/stats_${samplename}.txt

fi

if [ "$PAIRED" = "yes" ]
then

	if echo ${FILE} | grep -q "\.R1"
	then

			#calculate unique number of sequences
			echo -n "Unique number of sequences | " >> $OUTDIR/Quality/stats_${samplename}.txt
			dir=$( dirname $FILE )
			samplename_paired=$( get_sample_name ${FILE} | sed 's/\.R1//' )
			# paste of the two paired end fastq files
			paste <(zcat ${dir}/${samplename_paired}.R1.fastq.gz) <(zcat ${dir}/${samplename_paired}.R2.fastq.gz) > ${dir}/${samplename_paired}_paste.fastq
			# count the number of unique reads by paste the .R1 and .R2 reads and see if the concatenate sequence is unique
			grep -A 1 "@HWI" "${dir}/${samplename_paired}_paste.fastq" | sed "/@HWI/d" | sed "s/\t//g" | sed "/--/d" | sort | uniq | wc -l >> $OUTDIR/Quality/stats_${samplename}.txt
			rm "${dir}/${samplename_paired}_paste.fastq"

	elif echo ${FILE} | grep -q "\.R2"
	then
		
		samplename_paired=$( get_sample_name ${FILE} | sed 's/\.R1//' )
		#calculate unique number of sequences
		Unique_Number=$(grep "Unique number of sequences" $OUTDIR/Quality/stats_${samplename_paired}.R1.txt | cut -d"|" -f 2 | sed "s/ //")
		echo "Unique number of sequences | ${Unique_Number}" >> $OUTDIR/Quality/stats_${samplename}.txt
	fi

elif [ "$PAIRED" = "no" ]
then
	#calculate unique number of sequences
	echo -n "Unique number of sequences | " >> $OUTDIR/Quality/stats_${samplename}.txt
	zgrep -A 1 "@HWI" $FILE | sed "/@HWI/d" | sed "s/\t//g" | sed "/--/d" | sort | uniq | wc -l >> $OUTDIR/Quality/stats_${samplename}.txt
fi

#launch fastqc
echo "Launch FastQC..."
$BIN/$FASTQC_VERSION/fastqc --nogroup -o $OUTDIR/Quality/fastqc $FILE

wait 

echo "#----------------------------------------------------------"
echo "`date`: Check FASTQC version and unzip outputs if >= 0.11.2"
echo "#----------------------------------------------------------"

	
if [ ! -d "$OUTDIR/Quality/fastqc/${samplename}_fastqc/" ]
then

	if [ -f "$OUTDIR/Quality/fastqc/${samplename}_fastqc.zip" ]
	then

		unzip "$OUTDIR/Quality/fastqc/${samplename}_fastqc.zip" -d "$OUTDIR/Quality/fastqc/"

	else

		print_err_and_exit "File $OUTDIR/Quality/fastqc/${samplename}_fastqc.zip not found. Check you FastQC bin or version. Exit RNAseqPipeline."

	fi

fi


#draw 2 graphs :
#boxplot of quality score in each sequencing cycle
#% of each nt in each sequencing cycle
echo "Draw quality plots..."
fastqcfile=`echo $OUTDIR"/Quality/fastqc/"${samplename}"_fastqc/fastqc_data.txt"`

#line with sequence length information
seqlengthline=`grep -n "Sequence length" $fastqcfile | cut -f 1 -d ":"`

#first line with quality score information
qualityscoreline=`grep -n "Per base sequence quality" $fastqcfile | cut -f 1 -d ":"`

#first line with sequence content information
seqcontentline=`grep -n "Per base sequence content" $fastqcfile | cut -f 1 -d ":"`

#first line with N content information
Nseqcontentline=`grep -n "Per base N content" $fastqcfile | cut -f 1 -d ":"`
echo "source(\"$SCRIPTSDIR/src/Quality/DrawQualityPlotsFastQC0.9.5.R\");
filename = \"$fastqcfile\";title = \"${samplename}\";outdir = \"$OUTDIR/Quality\";
seqlength = $seqlengthline;qualityscore = $qualityscoreline;seqcontent = $seqcontentline;Nseqcontent = $Nseqcontentline;
qualityanalysis(filename, title, outdir, seqlength, qualityscore, seqcontent, Nseqcontent) " | $RBIN --no-save --no-restore --quiet
echo >> stats_${samplename}.txt
