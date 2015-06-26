#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
#This script performes the five parts of the auto report.
#
#If statistics option = no, in RNAseqPipeline, RNAseqReport performes
#these three steps:
#
#In a first part, it generate two tables which are include in the auto
#report:
#	- one which contain quality of sequencing experiment
#	- the other containing statistics of mapping
#
#In a second part, it generate a clustering graphe of all samples, which
#is include in the auto report.
#
#Finally, it generate a second report, named DataExploration, which
#contains different statitics on samples and generate different graphics
#which are include in the DataExploration report. Moreover, at this time,
#the Ensembl Version use for statistical analyses is checked. If it not 
#the same that use in RNAseq pipeline, the script crash and an error is
#reported.
#
#If statistics option = no, in RNAseqPipeline, RNAseqReport performes
#two supplementary step (others statistical analysis on the samples):
#The use of this script is optionnal in the RNAseq pipeline, because it 
#is necessary to know the different samples comparisons that we want to
#do in the experiment.
#
#In a first part, it generates a genes_diff.txt file of the experiment.
#This file contains the raw and normalized data and also the log2 Fold-
#Change and p-value for each comparisons.
#
#In a second part, it generates a scatterplot for each comparison, which
#is include in the auto report.
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh RNAseqReport.sh -c config_file -o outdir -d datadir -d design_file_report [-h]

		Description:
			This script performes an auto report on the RNAseq data.

		Options:
		 	-o, --outdir Outdir of the run
		        This option is required.
		        Generally of type : /path/to/Num_project_Name_author
		    -d, --datadir Datadir of the run
		    	This option is required. 
		    -c, --config-file Config file of the run.
		        This option is required. 
		    -f, --design-file-report Design file to use for the report.
		        This option is required.
		    -h, --help
		        Print this message and exit the program.
		__EOF__
}
#------------------------------------------------------------------------------------------------------------------

###################################################################################################################
# Getting parameters from the input
#------------------------------------------------------------------------------------------------------------------
ARGS=$(getopt -o "c:o:d:f:h" --long "config-file:,outdir:,datadir:,design_file_report:,help" -- "$@" 2> /dev/null)
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
		-d|--datadir)
			DATADIR="$2"
			shift 2 
			;;
		-f|--design_file_report)
			DESIGN_FILE_REPORT="$2"
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
[ "$DATADIR" == "" ] ||[ "$CONFIG_FILE" == "" ] || [ "$OUTDIR" == "" ] && \
	echo "Options --datadir, --config-file and --outdir are required. " "Use -h or --help to display the help." && exit 1;
#------------------------------------------------------------------------------------------------------------------

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

make_directories ${OUTDIR}/Report

cd ${OUTDIR}/Report

echo "`date` : Statistical analysis"

#Path R scripts (RNAseqAnalysis)
R_scripts=${SCRIPTSDIR}/src/StatisticalAnalysis/${SLURM_PARTITION}

#Copy of PRINCIPAL_FILE.tex: in this file be included all the .tex files necessary to generate auto report.
report_file="${OUTDIR}/Report/${PROJECT_NUMBER}_report.tex"

#The file containing quality table of the data
file_quality="${OUTDIR}/Quality/stats_summary.txt"

#The file containing the table of mapping quality of the data
file_mapping="${OUTDIR}/Alignment/tophat_stats_summary.txt"

DESIGN_FILE="$DESIGN_FILE_REPORT"

design_file="$OUTDIR/Scripts/design_file_modif.txt"
#######################################
# ADAPTING OF DESIGN FILE
#######################################  
#The file containing the design of the data analysis
cp ${DESIGN_FILE} ${design_file}
### view all lines except the first
		# tail -n +2 design_file.txt
### replace \t by any symbol
		# sed -e 's/\t/???/'
### replace the second \t by ",\t
		# sed -e 's/\t/",\t/'
### replace any symbol by \t"
		#sed -e 's/???/\t"/'
### you want to delete the last occurrence of a comma, for that we reverse the file
### it removes the first occurrence of the comma and the file is put back in order
		# sed -e '1!G;h;$!d;0,/,/ s///' | sed -e '1!G;h;$!d'
file_design="$OUTDIR/Scripts/design_fileOK.txt"
tail -n +2 ${design_file} | sed -e 's/\t/???/' | sed -e 's/\t/",\t/' | sed -e 's/???/\t"/' | sed -e '1!G;h;$!d;0,/,/ s///' | sed -e '1!G;h;$!d' > ${file_design}

head=$( head -1 ${design_file} )
echo ${head} > ${design_file}
cat ${file_design} >> ${design_file}

head -1 $design_file | sed -r 's/^ *//;s/ {1,}/\t/g' > ${file_design}
cat $design_file | sed 1,1d >> ${file_design}

cp ${file_design} ${design_file}
#######################################
#######################################  

#Recover of samplenames
samplename="${OUTDIR}/samplenames.txt"

#The file containing the design of the data analysis, use in R.
#The difference between design_file and design_file_R is that samplenames are separated by commas in the first 
#(in order to create the samplename in bash script) and the other not (to the analysis in R)

design_file_R="$OUTDIR/Scripts/design_file_R.txt"

#Recover of Bowtie version
bowtie=$( echo ${BOWTIE_VERSION} | grep -o "[^/]*$" | sed -e 's/_/\\_/g' )

#Recover of Tophat version
tophat=$( echo ${TOPHAT_VERSION} | grep -o "[^/]*$" | sed -e 's/.Linux.*//g' | sed -e 's/_/\\_/g' )

#Recover of Htseq version
htseq=$( echo ${HTSEQ_VERSION} | cut -d/ -f2  )

#The .tex file containing the part of quality and mapping analysis for auto report.
mapping_tex="${OUTDIR}/Report/${PROJECT_NUMBER}_Mapping.tex"

#The .tex file containing the part of clustering for auto report.
clustering_rnw="${OUTDIR}/Report/${PROJECT_NUMBER}_Clustering.rnw"
clustering_tex="${OUTDIR}/Report/${PROJECT_NUMBER}_Clustering.tex"

#The .tex file generating a supplementary report for data analysis
data_exploration_rnw="${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataExploration.rnw"
data_exploration_tex="${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataExploration.tex"

#Recovery of gene ensemble genome to use for annotation
genome_ensembl=$( grep -m 1 ${GENOME_VERSION} ${GENOME_CONVERSION} | cut -f2 )

#Path to report files of RNAseqDataStatistics
statistics_rnw="${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataStatistics.rnw"
statistics_tex="${OUTDIR}/Report/${PROJECT_NUMBER}_RNAseqDataStatistics.tex"

#Path to file genes_diff.txt
genes_diff=${OUTDIR}/Report/files/${PROJECT_NUMBER}_genes_diff.txt
#path to file genes_diff_DESeq2.txt
genes_diff_DESeq2=${OUTDIR}/Report/files/${PROJECT_NUMBER}_genes_diff_DESeq2.txt
allres_DESeq2=${OUTDIR}/Report/files/${PROJECT_NUMBER}_allres_DESeq2.txt
#Path to file allres.txt
allres_file=${OUTDIR}/Report/files/${PROJECT_NUMBER}_allres.txt
#Path to file sizeFactor.txt
sizeFactor=${OUTDIR}/Report/files/${PROJECT_NUMBER}_sizeFactor.txt

genome=$( echo $SPECIES | sed -e 's/\_/\\_/' )

eval GENE_BIOTYPE_FILE=$GENE_BIOTYPE_FILE

echo "### PARAMETERS FOR RNA-SEQ REPORT ###
OUTDIR : ${OUTDIR}
datadir : ${DATADIR}
Bowtie version : ${BOWTIE_VERSION}
Tophat version : ${TOPHAT_VERSION}
Genome version : ${GENOME_VERSION}
Gtf file : ${GTF_FILE}
Path to scripts of auto report : ${REPORTDIR}
Path to R : ${RBIN}
Project number : ${PROJECT_NUMBER}
Ensembl genome : ${genome_ensembl}
Statistical Analysis : ${STATISTICAL_ANALYSIS}
dist_method : ${DIST_METHOD}
hclust_method : ${HCLUST_METHOD}
SCRIPTSDIR : ${SCRIPTSDIR}
HTSEQ_VERSION : ${HTSEQ_VERSION}
species : ${SPECIES}
######
"

##################
# WITHOUT OPTION #
##################

###################
# START OF SCRIPT #
###################

touch ${report_file}
cp ${PRINCIPAL_FILE} ${report_file}

#The project number is written in the report
sed -i "s/PROJECT/${PROJECT_NUMBER}/" ${report_file}
#Insertion of the cover page
sed -i 's/%FrontCover/\\FrontCover/' ${report_file}

#####################
#     FUNCTIONS     #
#####################


#####################
# QUALITY & MAPPING #


########################################################################
#This function performes the first part RNAseqReport
#
#It generate two tables which are include in the auto
#report:
#
#	- one which contain quality of sequencing experiment
#	- the other containing statistics of mapping
#
########################################################################

########################################################################
#
# Parameters of this function:
#
# ${OUTDIR} ${samplename} ${file_quality} ${file_mapping} ${mapping_tex} 
# ${genome} ${GENOME_VERSION} ${tophat} ${bowtie} ${PROJECT_NUMBER}
#
########################################################################

QualityAndMapping(){

#Declaration of variables

OUTDIR=$1
samplename=$2
file_quality=$3
file_mapping=$4
mapping_tex=$5
genome=$6
GENOME_VERSION=$7
tophat=$8
bowtie=$9
PROJECT_NUMBER=${10}
PAIRED=${11}

assemblyOK=$( echo ${GENOME_VERSION} | sed -e 's/spikes//g' )

#Quality table modified for the report
quality_tab_tmp="${OUTDIR}/Quality/quality_tab.txt"
#Mapping table modified for the report
mapping_tab_tmp="${OUTDIR}/Quality/mapping_tab.txt"

#Concatenation of samplenames with tables and arrangement of the tables

if [ "$PAIRED" = "no" ]
then

	paste ${samplename} $file_quality > $quality_tab_tmp

# if PAIRED-end=yes, do not forget that there are two fastq files during qualityt analysis
# therefore must be doubled samplename
elif [ "$PAIRED" = "yes" ]
then
	samplename_tmp=${OUTDIR}/samplename_tmp.txt
	touch ${samplename_tmp}
	counter=0
	for ligne in `cat ${samplename}`
	do  
		if [ "$counter" -eq 0 ]
		then 
			echo $ligne >> ${samplename_tmp}
		else 
			echo $ligne >> ${samplename_tmp}
			echo $ligne >> ${samplename_tmp}
		fi
		
		counter=${counter+1}
	done  

	paste ${samplename_tmp} $file_quality > $quality_tab_tmp
	rm ${samplename_tmp}
	
fi
	
sed -i 's/"//g' $quality_tab_tmp
sed -i 's/File/Sample ID/' $quality_tab_tmp

paste ${samplename} $file_mapping > $mapping_tab_tmp
sed -i 's/"//g' $mapping_tab_tmp

#Creation of mapping.tex, which containing the two tables 
> ${mapping_tex}

#### PART OF QUALITY ANALYSIS ####

#Selection of columns.
#Tabs are replaced by "&" and line endings by "\\ \hline" to create the latex table.
Tabquality=$( cat $quality_tab_tmp | cut -f1,2,4,6,7,8 | sed 's/[\t]/ \& /g' | sed 's/$/ \\\\ \\hline/g' | sed 's/\_/\\_/g' | sed 's/\#/\\#/g' )

#Writting of quality table in mapping.tex file.
echo "
	
	\\vspace{1cm}
	\\textbf{\\Large{DATA}} \\newline
	A summary of the main characteristics of the data analyzed in this document is available in Table 1. \\newline
	\\begin{table}[H]
	\\hspace{-3cm}
	\\begin{tabular}{|p{3cm}|p{3cm}|p{1.5cm}|p{2cm}|p{3cm}|p{3cm}|} \\hline
	$Tabquality
	\\end{tabular}
	\\caption{Main characteristics of the data analyzed in this document}
	\\label{QUAL}
	\\end{table}
	\\vspace{0.25cm} 

" >> ${mapping_tex}

#### PART OF MAPPING  ####

#Selection of columns.
#Tabs are replaced by "&" and line endings by "\\ \hline" to create the latex table.
Tabmapping=$( cat $mapping_tab_tmp | cut -f1,2,3,6,7,8,9,10,11 | sed 's/[\t]/ \& /g' | sed 's/$/ \\\\ \\hline/g' | sed 's/\_/\\_/g' | sed 's/\#/\\#/g' | sed 's/\%/\\%/g' )

echo "
	\\newpage
	\\textbf{\\Large{MAPPING}} \\newline	
	Reads were mapped onto the ${assemblyOK} assembly of the ${genome} genome by using Tophat version : ${tophat} (\\cite{TRA2009}) and the bowtie version : ${bowtie} (\\cite{LAN2009}). A summary of the mapping results is available table in \ref{MAP}. Only uniquely aligned reads have been retained for further analyses. \\newline
	\\begin{table}[H]
	\\hspace{-3cm}
	\\begin{tabular}{|p{2.5cm}|p{2.5cm}|p{1.5cm}|p{2cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|p{1.5cm}|} \\hline
	$Tabmapping
	\\end{tabular} \\newline
	\\caption{Summary of mapping results}
	\\label{MAP}
	\\end{table}
	\\vspace{0.25cm}
	For each sample, an alignment file in BAM format (sampleid\\_alignment.bam) and the corresponding index (sampleid\\_alignment.bam.bai) are available on the following ftp server: ftp://ngs.igbmc.fr/analyzeddata.  \\newline
	\\newpage

" >> ${mapping_tex}

}

#####################
#     CLUSTERING    #

########################################################################
# This function performes the second part of RNAseqReport
# It generate a clustering graphe of all samples, which is include in 
# the auto report.
########################################################################

########################################################################
#
# Parameters of this function
#
# ${OUTDIR} ${design_file} ${clustering_rnw} ${clustering_tex}
#
########################################################################

#########################
# VARIABLES DECLARATION #
#########################


RNAseqDataClustering (){
	
#Declaration of variables	

OUTDIR=$1
design_file=$2
clustering_rnw=$3
clustering_tex=$4
PROJECT_NUMBER=$5
R_scripts=$6
RBIN=$7
HTSEQ_VERSION=$8

#Recovery of samplenames for use in R
samplenames=$( sed "1d" $design_file | cut -f2 )
echo ${samplenames} > ${OUTDIR}/samplename.txt
samplenames=$( cat ${OUTDIR}/samplename.txt )

###################
# START OF SCRIPT #
###################

#############
#START CHUNK#
#############
# Generates clustering graphic of all samples with the Pearson method.
# Creation of the clustering graphic is performed with the script 'RNAseqDataClustering.R'
# Other R commands, made ​​outside the script 'RNAseqDataClustering.R', allow to obtain the rawdata with samplenames.


echo "<<clustering, echo = FALSE, fig=TRUE>>=" >> $clustering_rnw	

echo "setwd(\"${OUTDIR}/Report\")" >> $clustering_rnw	

echo "samplename = c(${samplenames})" >> $clustering_rnw	

echo "datafile = read.table(\"${OUTDIR}/Quantification/htseq_combined.txt\", row.names=1, h=T, check.names=F, sep=\"\t\")" >> $clustering_rnw	

echo "rawdata = datafile[apply(datafile,1,sum)!=0,] " >> $clustering_rnw	

echo "colnames(rawdata) = samplename " >> $clustering_rnw	

echo "source(\"${R_scripts}/RNAseqDataClustering.R\")" >> $clustering_rnw	

echo "RNAseqClustering(rawdata, \"${DIST_METHOD}\", \"${HCLUST_METHOD}\")" >> $clustering_rnw	

echo "@" >> $clustering_rnw

###########
#END CHUNK#
###########

echo "\hspace{1cm} \\newline
	  \end{document}" >> $clustering_rnw

#Creation of clustering graphic with Sweave.

${RBIN} CMD Sweave ${clustering_rnw}

#Creation of pdf containing the clustering graphic.
${RBIN} CMD pdflatex ${clustering_tex}

#Creation of a new Clustering.tex file to insert the graphic and write comments.
> ${clustering_tex}

#Write of the results of clustering in Clustering.tex.
echo "
	\\vspace{1cm}
	\\textbf{\\Large{QUANTIFICATION}} \\newline

	Quantification of gene expression has been performed using ${htseq} (S.Anders \\footnote{http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html}). Read counts have then been normalized across libraries with the method proposed by Anders and Huber (\\cite{AND2010}).

	Following figure provides a clustering of all samples based on the correlation of gene expression between samples.

	\begin{figure}[!h]
	\centering
	\includegraphics[width=12cm]{${OUTDIR}/Report/${PROJECT_NUMBER}_Clustering-clustering.pdf}
	\caption{Clustering of all samples.
	Hierarchical clustering of all samples calculated using the UPGMA algorithm (Unweighted Pair Group Method with Arithmetic Mean) with the distance measure d=1-r (r=Pearson correlation coefficient).
	}
	\label{Clustering}
	\end{figure}
	
" > ${clustering_tex}

}

#####################
#    EXPLORATION    #

RNAseqDataExploration (){

########################################################################
# This function performes the third part of RNAseqReport
# It generates a supplementary report relative to auto report.
# This report contains different statistics on samples, for various
# checks of the analysis.
# Moreover, it contains different graphics on samples.
# The supplemantary repport is: RNAseqDataExploration.pdf
########################################################################

########################################################################
#
# Parameters of this function:
#
# ${OUTDIR} ${data_exploration_rnw} ${data_exploration_tex}
# ${R_scripts} ${RBIN} ${GTF_FILE} ${design_file} ${sizeFactor} ${PAIRED}
#
########################################################################

#Declaration of variables	

OUTDIR=$1
data_exploration_rnw=$2
data_exploration_tex=$3
R_scripts=$4
RBIN=$5
GTF_FILE=$6
design_file=$7
sizeFactor=$8
PAIRED=$9

#Recovery of samplenames for use in R
samplenames=$( cat "${OUTDIR}/samplename.txt" )

###################
# START OF SCRIPT #
###################

#############
#START CHUNK#
#############

gtfVersion=$( echo ${GTF_FILE} | grep -o "[^/]*$" )
host=$( grep ${gtfVersion} $HOST_CONVERSION | cut -f2 )

if [ -z ${host} ]
then
	echo " /!\ /!\ WARNING /!\ /!\ : the 'host_conversion.txt' file is not updated => Report failed."
fi

echo "\\newpage" >> ${data_exploration_rnw}
echo "\\textbf{\\Large{RNA-SEQ DATA EXPLORATION}} \\newline" >> ${data_exploration_rnw}

# Various graphics and statistics are performed with the script 'RNAseqDataExploration.R'
# Other R commands, made ​​outside the script 'RNAseqDataExploration.R'
# allow to obtain the rawdata with samplenames, ensemblVersion and load R libraries

echo "<<exploration, echo = FALSE>>=" >> ${data_exploration_rnw}	
	echo "setwd(\"${OUTDIR}/Report\")" >> ${data_exploration_rnw}
	echo "suppressMessages(library(DESeq))" >> ${data_exploration_rnw}
	echo "samplename = c(${samplenames})" >> ${data_exploration_rnw}
	echo "suppressMessages(library(biomaRt))" >> ${data_exploration_rnw}
	echo "martList = listMarts(host=\"${host}\")" >> ${data_exploration_rnw}
	echo "ensemblVersion = as.character(martList[martList['biomart']=='ENSEMBL_MART_ENSEMBL','version'])" >> ${data_exploration_rnw}
	echo "ensembl=useMart(\"ENSEMBL_MART_ENSEMBL\", host=\"${host}\", dataset='${genome_ensembl}')" >> ${data_exploration_rnw}
	echo "cat(ensemblVersion)" >> ${data_exploration_rnw}
	echo "datafile = read.table(\"${OUTDIR}/Quantification/htseq_combined.txt\", row.names=1, h=T, check.names=F, sep="	")" >> ${data_exploration_rnw}
	echo "rawdata = datafile[apply(datafile,1,sum)!=0,]" >> ${data_exploration_rnw}	
	echo "colnames(rawdata) = samplename" >> ${data_exploration_rnw}
	echo "nblib = dim(rawdata)[2]" >> ${data_exploration_rnw}
	# recovery of size factors to generate normalized wig files
  	echo "conds = factor(1:nblib)" >> ${data_exploration_rnw}
  	echo "cds = newCountDataSet(rawdata, conds)" >> ${data_exploration_rnw}
  	echo "cds = estimateSizeFactors(cds)" >> ${data_exploration_rnw}
  	echo "designFile=read.table(\"${design_file}\", sep= \"\t\", header=T,check.names = FALSE)" >> ${data_exploration_rnw}
  	echo "sizeFactorTable=as.data.frame(sizeFactors(cds))" >> ${data_exploration_rnw}
  	echo "write.table(sizeFactorTable, \"${sizeFactor}\", row.name=designFile[,1], sep=\"\t\")" >> ${data_exploration_rnw}
	# if the statistical analysis is not required, the annotation file is still generated
	# if the statistical analysis is required, the annotation file will be generated during this
	if [ "$STATISTICAL_ANALYSIS" = "no" ]
	then
		echo "source(\"${R_scripts}/RNAseqDataFormatingEnsVersion.R\")" >> ${data_exploration_rnw}
		echo "alldata = RNAseqDataFormatingEnsVersion(\"${OUTDIR}/Quantification/htseq_combined.txt\", \"${genome_ensembl}\", \"${host}\", samplename, \"${GENE_BIOTYPE_FILE}\", \"${GTF_FILE}\")" >> ${data_exploration_rnw}
		echo "write.table(alldata, \"${OUTDIR}/Report/files/${PROJECT_NUMBER}_raw_and_normalized_data.txt\", sep = \"\t\", row.names=FALSE)" >> ${data_exploration_rnw}
	fi
	echo "source(\"${R_scripts}/RNAseqDataExploration.R\")" >> ${data_exploration_rnw}
	echo "source(\"${SCRIPTSDIR}/src/Utilities/sere.R\")" >> ${data_exploration_rnw}
	echo "RNAseqDataExploration(rawdata, \"${sizeFactor}\", \"${design_file}\")" >> ${data_exploration_rnw}
	echo "sessionInfo()" >> ${data_exploration_rnw}
echo "@" >> ${data_exploration_rnw}

###########
#END CHUNK#
###########

echo "\end{document}" >> ${data_exploration_rnw}

cd ${OUTDIR}/Report

# We Start sweave on DataExploration file, a.tex is generated but also different graphics in folder "figures"
${RBIN} CMD Sweave ${data_exploration_rnw}

#Creation of pdf containing the various statistics.
${RBIN} CMD pdflatex ${data_exploration_tex}

sed -i 's/\\end{document}//' ${data_exploration_tex}

#Integration of graphics generated.

for i in `ls ${OUTDIR}/Report/figures/*.png`
do
	echo " \includegraphics[width=15cm]{${i}}" >> ${data_exploration_tex}
	
done

if [ "$PAIRED" = "yes" ]
then

	for i in `ls ${OUTDIR}/Quality/RSeQC/innerDistance/*.pdf`
	do
		echo " \includegraphics[width=15cm]{${i}}" >> ${data_exploration_tex}
		
	done
	
fi

wait

echo "\end{document}" >> ${data_exploration_tex}
	
} 

#####################
# QUALITY & MAPPING #
#####################

## grep du htseq et mise dans le bon ordre du design file

#sampleID=$( head -1 ${OUTDIR}/Quantification/htseq_combined.txt | sed "s/GeneID\t//" )

#head -1 $file_design > ${datadir}/tmp_design.txt

#for sample in $sampleID
#do
#	grep $sample $file_design >> ${datadir}/tmp_design.txt
#done

#mv ${datadir}/tmp_design.txt $file_design

#head -1 $design_file > ${datadir}/tmp_design.txt

#for sample in $sampleID
#do
#	grep $sample $design_file >> ${datadir}/tmp_design.txt
#done

#mv ${datadir}/tmp_design.txt $design_file



# Samplenames are recover from design file of the data analysis.
if [ -f $file_design -a -r $file_design ]
then
	
	#commas are removed to create the design file for use in R
	sed 's/,//' $design_file > $design_file_R
	#Take the column containing samplenames order to concatenate with the table of quality and mapping
	cut -f2 $design_file_R > $samplename
	#removing of empty lines in samplenames file
	sed -i '/^$/d' $samplename

else

	echo "Error with design file."
	exit 1
	
fi

#Check of files existence
if [ -f $file_quality -a -r $file_quality ]
then
    
	if [ -f $file_mapping -a -r $file_mapping ]
	then
		
		if [ -f $PRINCIPAL_FILE ]
		then
		
			#This function performes the first part of analysis, see the description above.
			#Parameters: ${OUTDIR}, ${samplename}, ${file_quality}, ${file_mapping}, ${mapping_tex}, ${genome}, ${GENOME_VERSION}
			# ${tophat}, ${bowtie}, ${PROJECT_NUMBER}
			
			QualityAndMapping ${OUTDIR} ${samplename} ${file_quality} ${file_mapping} ${mapping_tex} ${genome} ${GENOME_VERSION} ${tophat} ${bowtie} ${PROJECT_NUMBER} ${PAIRED}
			
			includeMapping="\\\input{${mapping_tex}}"
					
			#The mapping.tex file is included in the final report.
			
			sed -ie "s|%MAPPING|${includeMapping}|" ${report_file}
			
		else
			echo "Error with principal file of auto report."
			exit 1
		fi
	
	else
		echo "Error with file containing mapping statistics."
		exit 1
	fi

else
	echo "Error with file containing quality statistics."
	exit 1
fi

#####################
#    CLUSTERING     #
#####################

#This script Generates clustering graphic of all samples with the Pearson method.
#Parameters: ${OUTDIR} ${design_file} ${clustering_rnw} ${clustering_tex} ${PROJECT_NUMBER} ${R_scripts} ${RBIN}
		
#Creation of folders where results are stored
mkdir ${OUTDIR}/Report/figures ${OUTDIR}/Report/files

# ID of HTSeq_combined file are changed because R consider "Ensembl gene id" as 3 columns instead of one column.
sed -i 's/Ensembl gene id/GeneID/' ${OUTDIR}/Quantification/htseq_combined.txt
tr " " "\t" < ${OUTDIR}/Quantification/htseq_combined.txt > ${OUTDIR}/Quantification/htseq_combined_ok.txt
mv ${OUTDIR}/Quantification/htseq_combined_ok.txt ${OUTDIR}/Quantification/htseq_combined.txt

cd ${OUTDIR}/Report

cp ${PRINCIPAL_FILE} ${clustering_rnw}

sed -i 's/\\end{document}//' ${clustering_rnw}

RNAseqDataClustering ${OUTDIR} ${design_file} ${clustering_rnw} ${clustering_tex} ${PROJECT_NUMBER} ${R_scripts} ${RBIN} ${HTSEQ_VERSION}

#The Clustering.tex file is included in the final report.

includeClustering="\\\input{${clustering_tex}}"

sed -ie "s|%QUANTIFICATION|${includeClustering}|" ${report_file}

#####################
#  DATA EXPLORATION #
#####################

#The results are not included in the auto report but this script generates a supplementary report.
#This report contains different statistics on samples, for various checks of the analysis.
#Moreover, it contains different graphics on samples.
#Parameters: ${OUTDIR} ${data_exploration_rnw} ${data_exploration_tex} ${R_scripts} ${RBIN}

cp ${PRINCIPAL_FILE} ${data_exploration_rnw}

sed -i 's/\\end{document}//' ${data_exploration_rnw}

RNAseqDataExploration ${OUTDIR} ${data_exploration_rnw} ${data_exploration_tex} ${R_scripts} ${RBIN} ${GTF_FILE} ${design_file} ${sizeFactor} ${PAIRED}

#The project number is written in the report
sed -i "s/PROJECT/${PROJECT_NUMBER}/" ${data_exploration_tex}
#Insertion of the cover page
sed -i 's/%FrontCover/\\FrontCover/' ${data_exploration_tex}

${RBIN} CMD pdflatex ${data_exploration_tex}

#####################
#   END OF SCRIPT   #
#####################

mkdir ${OUTDIR}/Report/files/tex_files

rm *.texe

##################
#  WITH OPTION   #
##################

if [ "$STATISTICAL_ANALYSIS" = "yes" ]
then

	#########################
	#  STATISTICAL ANALYSIS #
	#########################

	#############
	#START CHUNK#
	#############

	cd ${OUTDIR}/Report

	cp ${PRINCIPAL_FILE} ${statistics_rnw}

	sed -i 's/\\end{document}//' ${statistics_rnw}

	#Write of the results of clustering in Clustering.tex.
	echo "

	\\newpage
	\\textbf{\\Large{STATISTICAL ANALYSIS}} \\newline
	Comparisons of interest have been performed using the statistical method proposed by Anders and Huber (\\cite{AND2010}). Resulting p-values were adjusted for multiple testing by using the Benjamini and Hochberg (\\cite{BEN1995}) method.
	The following figures (scatterplots) represents the comparisons of interest.
	The ${PROJECT_NUMBER}\_genes\_diff.txt file provides read counts and normalizes read counts for each gene together with the p-value and log2 Fold-change for the comparisons of interest.
	If we choose the thresholds False Discovery Rate (FDR) < 0.05 and |log2 Fold-change| > 1, the following number of significantly differentially expressed genes are found in each comparison: \footnote{In a comparison between A and B, an overexpressed gene is a gene with an higher expression in condition B than in condition A (log2 (B/A) > 1.} \newline
		
	" >> ${statistics_rnw}

	#########################################################################################################################
	# R commands, made ​​outside the scripts, allow to obtain the rawdata with samplenames, 
	# ensembl version and load R libraries or R scripts
	# The script 'RNAseqDataFormating.R' performed 'alldata' from htseq_combined.txt
	# alldata contains annotation of all genes and raw and normalized reads counts for each gene
	# The script 'RNAseqDataStatistics.R' performes all the comparisons described in design file
	# It uses "comp2cond.R" script if there are replicates in the experiment
	# and it uses "comp2cond_var2samples.R" script if there are no replicates
	# The script 'RNAseqDataStatistics.R' generates numProject_genes_diff.txt, which contains the same columns that alldata
	# but with supplementary columns: the log2 fold-change & p-value for each comparison
	# Moreover, it generates allres.txt, which allow to create venn diagramms 
	# with overexpressed and underexpressed genes between two samples
	# Finally, this script generate one scatterplot per comparison, with deregulated genes in red
	#########################################################################################################################


	 #For other genome version (e.g. version 67 for mouse): martList = listMarts(host=\"may2012.archive.ensembl.org/biomart/martservice/\")
	 #ensemblVersion = as.character(martList[martList['biomart']=='ENSEMBL_MART_ENSEMBL','version'])
	 #ensembl=useMart(\"ENSEMBL_MART_ENSEMBL\", host=\"may2012.archive.ensembl.org/biomart/martservice/\", dataset='mmusculus_gene_ensembl')
	 
	 # recover of genome version to use
	gtfVersion=$( echo ${GTF_FILE} | grep -o "[^/]*$" )
	host=$( grep ${gtfVersion} $HOST_CONVERSION | cut -f2 )

	echo "<<stats, echo = FALSE, fig=FALSE>>=" >> ${statistics_rnw}

	echo "setwd(\"${OUTDIR}/Report\")" >> ${statistics_rnw}
	echo "suppressMessages(library(DESeq))" >> ${statistics_rnw}
	echo "samplename = c(${samplenames})" >> ${statistics_rnw}
	echo "suppressMessages(library(biomaRt))" >> ${statistics_rnw}
	echo "martList = listMarts(host=\"${host}\")" >> ${statistics_rnw}
	echo "ensemblVersion = as.character(martList[martList['biomart']=='ENSEMBL_MART_ENSEMBL','version'])" >> ${statistics_rnw}
	echo "ensembl=useMart(\"ENSEMBL_MART_ENSEMBL\", host=\"${host}\", dataset='${genome_ensembl}')" >> ${statistics_rnw}
	echo "source(\"${R_scripts}/RNAseqDataFormatingEnsVersion.R\")" >> ${statistics_rnw}
	echo "alldata = RNAseqDataFormatingEnsVersion(\"${OUTDIR}/Quantification/htseq_combined.txt\", \"${genome_ensembl}\", \"${host}\", samplename, \"${GENE_BIOTYPE_FILE}\", \"${GTF_FILE}\")" >> ${statistics_rnw}
	echo "datafile = read.table(\"${OUTDIR}/Quantification/htseq_combined.txt\", row.names=1, h=T, check.names=F, sep="	")" >> ${statistics_rnw}
	echo "rawdata = datafile[apply(datafile,1,sum)!=0,]" >> ${statistics_rnw}	
	echo "colnames(rawdata) = samplename" >> ${statistics_rnw}
	echo "write.table(alldata, \"${OUTDIR}/Report/files/${PROJECT_NUMBER}_raw_and_normalized_data.txt\", sep = \"\t\", row.names=FALSE)" >> ${statistics_rnw}
	echo "source(\"${R_scripts}/comp2cond.R\")
	source(\"${R_scripts}/comp2cond_var2samples.R\")" >> ${statistics_rnw}
	echo "source(\"${R_scripts}/comp2cond_DESeq2.R\")" >> ${statistics_rnw}

	echo "source(\"${R_scripts}/RNAseqDataStatistics.R\")
	source(\"${R_scripts}/RNAseqDataStatistics_DESeq2.R\")" >> ${statistics_rnw}

	echo "RNAseqDataStatistics(\"${design_file}\", rawdata, \"${genes_diff}\", \"${allres_file}\", \"${THRESHOLD_LOG_FC}\", \"${THRESHOLD_ADJ_PVAL}\")" >> ${statistics_rnw}


	if [ "${SLURM_PARTITION}" = "rufus" ]
	then

		echo "RNAseqDataStatisticsDESeq2(\"${design_file}\", rawdata, alldata, \"${genes_diff_DESeq2}\", \"${THRESHOLD_LOG_FC}\", \"${THRESHOLD_ADJ_PVAL}\", \"${FIT_TYPE}\")" >> ${statistics_rnw}

	else

		echo "RNAseqDataStatisticsDESeq2(\"${design_file}\", rawdata, alldata, \"${genes_diff_DESeq2}\", \"${allres_DESeq2}\", \"${THRESHOLD_LOG_FC}\", \"${THRESHOLD_ADJ_PVAL}\", \"${FIT_TYPE}\")" >> ${statistics_rnw}

	fi

	echo "@" >> ${statistics_rnw}

	echo "\end{document}" >> ${statistics_rnw}
	###########
	#END CHUNK#
	###########

	#Creation of clustering graphic with Sweave.

	${RBIN} CMD Sweave ${statistics_rnw}

	wait

	#Creation of pdf containing the clustering graphic.
	${RBIN} CMD pdflatex ${statistics_tex}

	#Copy the contents of file S12050_RNAseqDataStatistics.tex in a variable because after sweave and pdflatex commands
	#we obtain a .tex containing \begin{document}, \end{document} etc ...
	#and they should be deleted for include S12050_RNAseqDataStatistics.tex in S12050_report.tex

	string=`cat ${statistics_tex} `

	#Creation of a new S12050_RNAseqDataStatistics.tex file to insert only the results of sweave
	#and not the structure of latex document necessary for sweave before (usepackage, \begin{document}...)
	> ${statistics_tex}

	#Removes all that is before REFERENCES in S12050_RNAseqDataStatistics.tex (corresponding to usepackages)
	echo "${string##*REFERENCES} " > ${statistics_tex}


	#Removes everything that begins with '\end' or '\begin'
	sed -i '/\\end/d' ${statistics_tex}
	sed -i '/\\begin/d' ${statistics_tex}

	sed -i 's/\\_/_/g' ${statistics_tex}
	sed -i 's/_/\\_/g'  ${statistics_tex}
	sed -i 's/\[1\]//g' ${statistics_tex}
	##############################
	# GENERATION OF SCATTERPLOTS #
	##############################

	#Create a folder named 'scatterplots' to store all scatterplots.png
	mkdir ${OUTDIR}/Report/figures/scatterplots
	mv ${OUTDIR}/Report/figures/*_vs_*.png ${OUTDIR}/Report/figures/scatterplots

	cd ${OUTDIR}/Report/figures/scatterplots/

	sed -i '/ERCC/d' "${OUTDIR}/Report/files/${PROJECT_NUMBER}_genes_diff.txt"
	sed "s/Inf/100000/g" "${OUTDIR}/Report/files/${PROJECT_NUMBER}_genes_diff.txt" > "${OUTDIR}/Report/files/${PROJECT_NUMBER}_genes_diff_without_inf.txt"
	sed -i -e "s/-Inf/-100000/g" "${OUTDIR}/Report/files/${PROJECT_NUMBER}_genes_diff_without_inf.txt"

	#A figure environment is created for each scatterplot finding in folder 'scatterplots'
	#Scatterplots will be include in auto report
	for i in `ls *scatterplot.png`
	do

		name=$( echo $i | cut -d. -f1 | sed 's/_vs_/ and /g' | sed 's/_/ /g' )

		echo "\begin{figure}[htbp]
			\centering
			\includegraphics[width=12cm]{${OUTDIR}/Report/figures/scatterplots/${i}}
			\caption{Comparison between ${name} (with DESeq).}
			\label{${name}}	\end{figure}" >> ${statistics_tex}
	done

	#Integration of each MA plot

	for i in `ls *_MA.png`
	do

		name=$( echo $i | cut -d. -f1 | sed 's/_vs_/ and /g' | sed 's/_/ /g' | sed 's/MA//g' )

		echo "\begin{figure}[htbp]
			\centering
			\includegraphics[width=12cm]{${OUTDIR}/Report/figures/scatterplots/${i}}
			\caption{MA plot between ${name} (with DESeq).}
			\label{${name}}	\end{figure}" >> ${statistics_tex}
	done

	for i in `ls MA_*.png`
	do

		name=$( echo $i | cut -d. -f1 | sed 's/_vs_/ and /g' | sed 's/MA_Condition_//g' | sed 's/_DESeq2//g' | sed 's/_/ /g'  )

		echo "\begin{figure}[htbp]
			\centering
			\includegraphics[width=12cm]{${OUTDIR}/Report/figures/scatterplots/${i}}
			\caption{MA plot between ${name} (with DESeq2).}
			\label{${name}  DESEQ2}	\end{figure}" >> ${statistics_tex}
	done

	#RNAseqDataStatistics.tex is included in auto report
	#It contains description of the analyses, statistics on differential expression and scatterplots
	includeStatistics="\\\input{${statistics_tex}}"

	sed -ie "s|%STATISTICS|${includeStatistics}|" ${report_file}

	sed -i 's/\\end{document}//' ${report_file}

	echo "
	\newpage
	\bibliographystyle{unsrt}
	\bibliography{${BIBTEX_FILE}}
	\end{document} " >> ${report_file}	
		
	sed -i 's/\\end{document}//' ${data_exploration_tex}

	echo	"\includegraphics[width=15cm]{${OUTDIR}/Report/figures/variance.pdf}
	\newline \centering{Variance between all samples.} \\" >> ${data_exploration_tex}
		
	cd ${OUTDIR}/Report/figures/
		
	for i in `ls *histpval.pdf`
	do

		name=$( echo $i | cut -d. -f1 | sed 's/_/ /g' | sed 's/histpval//g' )

		echo "\includegraphics[width=15cm]{${OUTDIR}/Report/figures/${i}} 
		\newline \centering{Comparison between ${name}} \\" >> ${data_exploration_tex}
		
	done
	echo "
	\newpage
	\bibliographystyle{plain}
	\nocite{*}
	\bibliography{${BIBTEX_FILE}} 
	\end{document}" >> ${data_exploration_tex}
	
	rm $file_design $design_file_R $design_file
fi

wait

rm ${OUTDIR}/samplenames.txt ${OUTDIR}/samplename.txt