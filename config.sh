#! /bin/bash
#------------------------------------------------------------------------------------------------------------------
# Authors: Amandine VELT
# Mail: velt@igbmc.fr
# Date: Avril 2015
#------------------------------------------------------------------------------------------------------------------
# This script contains variables that are used by the RNA-seq pipeline.
# Please check all these variables before launching the RNA-seq pipeline.
# It must be loaded with a source command.
###################################################################################################################

######################################################
# /!!!\ if you are on furious, don't forget to change $bin and $rbin params /!!!\ 
######################################################

###################################################################################################################
# Slurm variables
#------------------------------------------------------------------------------------------------------------------
# This pipeline is developed in order to use slurm workload manager to control the jobs.
# If the system does not use slurm please set the variable below to false. In this case the pipeline will be
# executed as a list of commands and will work slowly.
USE_SLURM=true

# Choose furious or rufus
SLURM_PARTITION="rufus"

# Command that list ids of all active jobs
SQUEUE="squeue | tr -s ' ' | cut -d ' ' -f 2"

#------------------------------------------------------------------------------------------------------------------
# To execute the command, the pipeline uses slurm sbatch command. Please give the absolute path of this software,
# otherwise it is assuming that sbatch is somewhere in your PATH. 
SBATCH=sbatch
SRUN=srun
# The following options are used when using sbatch. If an option is empty then it will be omitted if it is called.
# If you use another command than sbatch, please retrieve and set the equivalent options. Slurm definition of these 
# options are given to help you finding an equivalent.
#
# SLURM: "Specify a name for the job allocation. The specified name will appear along with the job id number when 
# querying running jobs on the system. The default is the name of the batch script, or just "sbatch" if the script 
# is read on sbatch’s standard input."
SBATCH_JOB_NAME="--job-name"
# SLURM: "Request a specific partition for the resource allocation. If not specified, the default behavior is to 
# allow the slurm controller to select the default partition as designated by the system administrator. If the 
# job can use more than one partition, specify their names in a comma separate list."
SBATCH_PARTITION="--partition"
# SLURM: "Instruct SLURM to connect the batch script’s standard output directly to the file name specified."
SBATCH_STDOUT="--output"
# SLURM: "Instruct SLURM to connect the batch script’s standard error directly to the file name specified."
SBATCH_STDERR="--error"
# SLURM: "Defer  the  start  of  this job until the specified dependencies have been satisfied completed.
# <dependency_list> is of the form <type:job_id[:job_id][,type:job_id[:job_id]]>.   Many  jobs  can share the same 
# dependency and these jobs may even belong to different   users.  The   value may be changed after job submission
# using the scontrol command.
# after:job_id[:jobid...]
#        This job can begin execution  after  the  specified  jobs
#        have begun execution.
# afterany:job_id[:jobid...]
#        This  job  can  begin  execution after the specified jobs
#        have terminated.
# afternotok:job_id[:jobid...]
#        This job can begin execution  after  the  specified  jobs
#        have terminated in some failed state (non-zero exit code,
#        node failure, timed out, etc).
# afterok:job_id[:jobid...]
#        This job can begin execution  after  the  specified  jobs
#        have  successfully  executed  (ran  to completion with an
#        exit code of zero).
# expand:job_id
#        Resources allocated to this job should be used to  expand
#        the specified job.  The job to expand must share the same
#        QOS (Quality of Service) and partition.  Gang  scheduling
#        of resources in the partition is also not supported.
# singleton
#        This   job  can  begin  execution  after  any  previously
#        launched jobs sharing the same job  name  and  user  have
#        terminated."
SBATCH_DEPENDENCY="--dependency"
SBATCH_DEPENDENCY_OPTION="afterok:"
SBATCH_DEPENDENCY_SEPARATOR=":"
# SLURM: "Request that ntasks be invoked on each node."
SBATCH_TASKS_PER_NODE="--ntasks-per-node"
# SLURM: "Set the working directory of the batch script to directory before it is executed. The path can be 
# specified as full path or relative path to the directory where the command is executed."
SBATCH_WORKING_DIRECTORY="--workdir"
# SLURM: "Specify the real memory required per node in MegaBytes."
SBATCH_MEMORY="--mem"
###################################################################################################################



###################################################################################################################
#information relative to the project
#------------------------------------------------------------------------------------------------------------------
#path/to/output directory where to store all the results, without final /
OUTDIR_ANALYSIS=$( echo "`date "+%Y-%m-%d %H:%M:%S"`" | sed "s/ /_/" ) #*#*#*
# name of data folder
DATA_FOLDER="rawdata"
#ID of the project
PROJECT_NUMBER="S14096" #*#*#*
###################################################################################################################



###################################################################################################################
#information on library preparation and sequencing
#------------------------------------------------------------------------------------------------------------------
#whether the library preparation was strand-specific. Values=no or reverse
STRANDED="reverse" #*#*#*
#whether spikes have been added during library preparation. Values=yes or no
SPIKES="yes" #*#*#*
#whether paired-end sequencing was performed. Values=yes or no
PAIRED="no" #*#*#*
# This is the expected (mean) inner distance between mate pairs => only use if paired=yes
FRAGMENT_SIZE="50"
###################################################################################################################



###################################################################################################################
#information on species, genome and annotations
#------------------------------------------------------------------------------------------------------------------
#name of the species used (eg Homo_sapiens or Mus_musculus) as specified in /ifs/illumina/share/Genomes directory 
SPECIES="Homo_sapiens" #*#*#*
#assembly version (eg mm10 or hg19) if there are spikes specify "spikes" after the assembly version (eg mm10spikes or hg19spikes)
GENOME_VERSION="hg19" #*#*#*
## path where all the genome's files are stocked
GENOME_PATH="/ifs/illumina/share/Genomes"
## path where all the genome's files are stocked
SPIKES_PATH="/ifs/illumina/share/Spikes"
#path to bowtie indexes
BOWTIE_INDEXES="${GENOME_PATH}/${SPECIES}/${GENOME_VERSION}/Bowtie"
#path to GTF file containing gene annotations
GTF_FILE="${GENOME_PATH}/${SPECIES}/${GENOME_VERSION}/Annotations/Ensembl/Homo_sapiens.GRCh37.69_UCSConlychr.gtf" #*#*#*
#path to bowtie index for spikes alignment
INDEX_SPIKES="${SPIKES_PATH}/Bowtie/spikes_only"
###################################################################################################################



###################################################################################################################
#options for the analyses performed by the pipeline
###################################################
#number of processors on which to launch tophat/tophat2 and some other tools
NBPROC="2" #*#*#*
#If yes, TopHat will map every read in all the mapping steps, reporting the best possible alignment found in any of these mapping steps. 
#This may greatly increase the mapping accuracy at the expense of an increase in running time.
# option only for tophat2 !!!
REALIGN="no"
#whether to perform statistical analysis or not. Values=yes or no
STATISTICAL_ANALYSIS="yes" #*#*#*

#options for statistical analyses
DIST_METHOD="pearson" #distance method to use for data clustering 
HCLUST_METHOD="average" #clustering method to use for data clustering
THRESHOLD_LOG_FC=1
THRESHOLD_ADJ_PVAL=0.05
FIT_TYPE="parametric"  #either "parametric", "local", or "mean" for the type of fitting of dispersions to the mean intensity.
#design of the analysis
DESIGN_FILE='$OUTDIR/Scripts/design_file.txt' # create via the ${SCRIPTSDIR}/RNAseqReport/Files/design_file.txt model file. This file containing the design of the analysis ( replicates, comparisons ... )
#window size to create the wig & tdf files => tipically "1" for RNA-seq
WINDOW_SIZE="1"
###################################################################################################################



###################################################################################################################
#scripts and software versions and location
###########################################
 #directory to store temporary files (caution, must have a lot of available disk space)
TMP_DIR="/ifs/home/velt/tmp/"
#path/to/bin directory. This prefix is used for all the following binary locations
BIN="/ifs/illumina/share/software" # for furious : /ifs/illumina/share/Utilities/softwareSL
#path/to/R binary for furious ${BIN}/R/R-3.1.1/bin/R
RBIN="/biolo/R_surf/R-3.0.1/bin/R" # for furious : ${BIN}/R/R-3.1.1/bin/R
#path/to/directory where to find all the scripts needed by this pipeline
SCRIPTSDIR="/ifs/home/velt/Git/RNAseqPipeline" #*#*#*
#path to Utilities scripts
UTILITIESDIR="/ifs/illumina/share/Utilities"
# path/to/directory where to find all the scripts needed to analyse spikes
SCRIPTSDIR_SPIKES="${SCRIPTSDIR}/src/Spikes"
#path/to/fastqc binary directory, without final / (i.e. fastqc binary must be located in : ${BIN}/fastqcbindir/fastqc)
FASTQC_VERSION="FastQC/fastqc_v0.11.2"
#path/to/tophat binary directory, without final / (i.e. tophat binary must be located in : ${BIN}/tophatbindir/tophat)
TOPHAT_VERSION="tophat/tophat-2.0.14.Linux_x86_64"
#path/to/tophat binary directory, without final / (i.e. tophat binary must be located in : ${BIN}/tophatbindir/tophat) 
BOWTIE_VERSION="bowtie/bowtie-2-2.1.0"
#path/to/samtools binary directory, without final / (i.e. samtools binary must be located in : ${BIN}/samtoolsbindir/samtools)
SAMTOOLS_VERSION="Samtools/samtools-0.1.19"
#path/to/htseq-count binary directory, without final / (i.e. htseq-count binary must be located in : ${BIN}/htseqbindir/htseq-count)
HTSEQ_VERSION="HTSeq/HTSeq-0.6.1/build/scripts-2.7"
#path/to/htseq-count python library => the .egg use by python to run htseq count
HTSEQ_PYTHONPATH="HTSeq-0.6.1-py2.7-linux-x86_64.egg"
#path/to/RSeQC binary directory, without final
RSEQC_VERSION="RSeQC/RSeQC-2.5/scripts" 
#path/to/RSeQC binary directory, without final
RSEQC_PYTHONPATH="RSeQC-2.5-py2.7-linux-x86_64.egg" 
#path/to/IGV binary directory, without final /
IGV_VERSION="IGV/IGVTools_2.2.1"
#path/to/directory where to find all the scripts needed to generate automatic report
REPORTDIR="${SCRIPTSDIR}/src/RNAseqReport"
# path to python 2.7 libraries
# pour htseq, on ajoute le nom de dossier du package python après ce pythonpath
# donc faire attention que l'appelation des dossier .egg d'htseq ne change pas
PYTHONPATH="${BIN}/Python/lib/python2.7/site-packages"
###################################################################################################################



###################################################################################################################
#data files used by the pipeline
################################
# principal .tex file for generation of .pdf report for spikes.
PRINCIPAL_FILE_SPIKES="${SCRIPTSDIR_SPIKES}/Report/principal_file.tex" 
# file containing the sequence length of each spike.
LENGTH_SPIKES="${SCRIPTSDIR_SPIKES}/Files/Spikes_Length_sort.txt" 
# file containing the concentration of each spike in each mix.
CONCENTRATION_SPIKES="${SCRIPTSDIR_SPIKES}/Files/Concentration.txt"
# file containing the expected fold-change for each spike between the two mix.
FOLDCHANGE_SPIKES="${SCRIPTSDIR_SPIKES}/Files/Expected_FoldChange.txt" 
SCRIPT_SPIKES_POSITIONS="${SCRIPTSDIR_SPIKES}/Positions_Alignment.pl"
# file containing information on all spike sequences
GTF_SPIKES="${SPIKES_PATH}/Annotations/ERCC92.gtf"
SPIKES_ANNOTATIONS="${SPIKES_PATH}/Annotations/ERCC92.bed"
# file containing the conversion of genome assembly version to ensembl genome name and the corresponding fasta file.
GENOME_CONVERSION="${SCRIPTSDIR}/src/RNAseqReport/Files/genome_conversion.txt" 
# file using gtf file to retrieve the ensembl version to use for annotations.
HOST_CONVERSION="${SCRIPTSDIR}/src/RNAseqReport/Files/host_conversion.txt" 
# principal .tex file for generation of .pdf report.
PRINCIPAL_FILE="${SCRIPTSDIR}/src/RNAseqReport/Files/principal_file.tex" 
# file containing the bibliography of the automated report.
BIBTEX_FILE="${SCRIPTSDIR}/src/RNAseqReport/Files/bibliography"
# file containing all the statistics results aubout all samples
STATISTICS_FILE='${OUTDIR}/Quality/all_statistics_results.txt'
# file containing statistics on gene biotype
GENE_BIOTYPE_FILE='${OUTDIR}/Report/files/geneBiotypeStats.txt'
###################################################################################################################



###################################################################################################################
#log files used by the pipeline
################################
# log file of automatic report for spikes
LOGFILE_SPIKES='${OUTDIR}/Logs/log_file_Spikes.txt'
# log file of automatic report + data exploration report
LOGFILE_REPORT='${OUTDIR}/Logs/log_file_Report.txt'
# log file of wig and tdf files generation
LOGFILE_WIG_TDF='${OUTDIR}/Logs/log_file_wigTdf.txt'
###################################################################################################################
