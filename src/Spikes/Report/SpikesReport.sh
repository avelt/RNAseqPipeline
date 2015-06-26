#! /bin/sh
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: March 2015
#------------------------------------------------------------------------------------------------------------------
#
# Parameters of this script: given by RNAseqPipelineAllsamples.sh
#
# ${outdir}, ${datadir}, ${projectNumber}, ${rbin}, ${designfile} 
# => recover via the parameter file
#
# Generation of an automatic report for spikes-in 
#
#------------------------------------------------------------------------------------------------------------------


###################################################################################################################
# Define the function to print the usage of the script
#------------------------------------------------------------------------------------------------------------------
function usage
{
	cat <<-__EOF__
		Usage:
		    sh RNAseqPipeline1sample.sh -c config_file -o outdir [-h]

		Description:
			Generation of an automatic report for spikes-in 

		Options:
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

source $CONFIG_FILE

source "$SCRIPTSDIR/src/Utilities/utils.sh"

# Copy of principal_file.tex: in this file be included the frame of auto report.
report_file="${OUTDIR}/Report/Spikes/${PROJECT_NUMBER}_spikes_report.tex"

#The file containing the table of mapping quality of the data
file_mapping="${OUTDIR}/Alignment/tophat_stats_summary.txt"

#Recover of samplenames
samplename="${OUTDIR}/samplenames.txt"

# temporary design file suitable for use in R
DESIGN_FILE="$OUTDIR/Scripts/design_file.txt"
design_file_R="$OUTDIR/Scripts/design_file.txt"

###################
# START OF SCRIPT #
###################

touch ${report_file}
cp ${PRINCIPAL_FILE_SPIKES} ${report_file}

#The project number is written in the report
sed -i "s/PROJECT/${PROJECT_NUMBER}/" ${report_file}

# Samplenames are recover from design file of the data analysis.
if [ -f $DESIGN_FILE -a -r $DESIGN_FILE ]
then
	
	#commas are removed to create the design file for use in R
	sed 's/,//' $DESIGN_FILE > $design_file_R
	#Take the column containing samplenames order to concatenate with the table of quality and mapping
	cut -f2 $design_file_R > $samplename
	#removing of empty lines in samplenames file
	sed -i '/^$/d' $samplename

else

	print_err_and_exit "Error with design file."
	
fi

#############
#FIRST STEP : Number of reads corresponding to spikes and  % of spikes represented by reads
#############

# Sort by samples ID to have the same order as in tophat table.
sort -k 1 ${OUTDIR}/Quality/Spikes/Number_of_reads.txt | cut -d: -f2 > $OUTDIR/tmp_spikes.txt

# Removes unnecessary spaces
sed -i 's/%/\%/g' $OUTDIR/tmp_spikes.txt

# Recovery of columns
TabSpikes=$( cat ${OUTDIR}/Quality/Spikes/Number_of_reads.txt | sed 's/\_/\\_/g' | sed 's/\%/\\%/g' | sed 's/\:/\&/g' | sed 's/$/ \\\\ \\hline/g' )

# Writing the array and the percent of reads representing spikes in "report_file"

echo "
	\\centering{\\textbf{\\Large{MAPPING}}} \newline
	
	A table showing the results of mapping. The last column represents the number of reads corresponding to spikes. Ideally, this number should be 1\% of total reads.
" >> ${report_file}

echo "
	\\begin{table}[H]
	\centering
	\\hspace{-3cm}
	\\begin{tabular}{|p{6.5cm}|p{6.5cm}|} \\hline
	${TabSpikes}
	\\end{tabular} \\newline
	\\caption{Summary of mapping results}
	\\label{MAP}
	\\end{table}
	\\textbf{$percent_spikes}

" >> ${report_file}

echo "
	(Note: Some spikes are in very low concentration, it is normal that 100\% of spikes are not represented by reads. However, must pay attention to what the majority is represented.)
\\newpage" >> ${report_file}

##############
#SECOND STEP : correlation coefficient between concentration and NRN
##############

cd ${OUTDIR}/Quality/Spikes

echo "
	\\centering{\\textbf{\\Large{CORRELATION COEFFICIENT}}} \newline
	
	Graphs showing the concentration of spikes according to NRN. 
	Each spike has a different concentration in the Spikes mix used. And here we want to know if the number of reads obtained for each spike is proportional to the concentration of this spike. In these graphs, there are 92 points, corresponding to 92 spikes.
	NRN : After normalization of reads number, it is divided by the size of the spike which it correspond to obtain the final number of reads normalized (NRN: Number of Reads Normalized) for each spike in each sample.
" >> ${report_file}

# Recovery of all NRN graphics (one per sample) and add in report_file
for i in `ls ${OUTDIR}/Quality/Spikes/*NRN.png`
do

#Some names have points, these points are suppressed and replaced by a dash. 
#Then renames the graph with the mv command. Otherwise latex generates error.
name1=$( echo $i | sed 's/_NRN.png//g' | sed 's/\./-/g' ) 
samplename=$( basename $name1 )
mv ${i} ${OUTDIR}/Quality/Spikes/${samplename}_NRN.png
#Here we replaced "_" by "\\_" for that it is compatible in latex, otherwise error
samplename_title=$( basename $name1 | sed 's/_/\\_/g' | sed 's/\./-/g' )



echo "\begin{figure}[H]
	\centering
	\includegraphics[width=12cm]{${OUTDIR}/Quality/Spikes/${samplename}_NRN.png}
	\caption{Coefficient correlation for ${samplename_title} sample .}
	\label{NRN ${name1}}	\end{figure} " >> ${report_file}
done

echo "\\newpage" >> ${report_file}

##############
#THIRST STEP : Fold-change 
##############

# Recovery of fold-change graphic and add in report_file

echo "
	\\centering{\\textbf{\\Large{FOLD-CHANGE}}} \newline
	
	Graph showing the expected fold-change according to observed fold-change. 
	The spikes kit is composed of two mixes of spikes, composed of the same spikes but in different concentration.
	The kit gives the values of fold-change expected between the two mixes for each spike. And this graph shows if the fold-change are well preserved during the experiment.
	Note: In this graph, there are 92 points, corresponding to 92 spikes and the red dots correspond to spikes represented by less than 50 reads.
" >> ${report_file}

echo "\begin{figure}[H]
	\centering
	\includegraphics[width=12cm]{${OUTDIR}/Quality/Spikes/foldChange_spikes.png}
	\caption{Comparison between expected and observed fold-change for mix 1 and 2 of spikes-in.}
	\label{fold-change}	\end{figure}" >> ${report_file}

echo "\\newpage" >> ${report_file}
############
#FOUR STEP : Coverage
############

# Recovery of all coverage graphics (one per sample) and add in report_file

echo "
	\\centering{\\textbf{\\Large{COVERAGE ALONG SPIKES}}} \newline
	
	Graphs showing number of reads at each position of spikes (size of the spikes is represented in percentage).
	This allows to see if the coverage of reads along the spikes is homogeneous or not, if there are biases ...
" >> ${report_file}

if [ -f "${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples.geneBodyCoverage.curves.pdf" ]
then

mv ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples.geneBodyCoverage.curves.pdf ${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples_geneBodyCoverage_curves.pdf

png="${OUTDIR}/Quality/Spikes/geneBodyCoverage/allsamples_geneBodyCoverage_curves.pdf"
name1=$( basename $png )

echo "\begin{figure}[H]
	\centering
	\includegraphics[width=12cm]{${png}}
	\caption{Coverage of reads along spikes for all samples.}
	\label{Coverage ${name1}}	\end{figure}" >> ${report_file}

else

	echo "Gene body coverage file doesn't exist so doesn't include in the Spikes report !!"

fi

echo "\end{document}" >> ${report_file}

rm ${OUTDIR}/tmp_spikes.txt

cd ${OUTDIR}/Report/Spikes

${RBIN} CMD pdflatex ${report_file}
