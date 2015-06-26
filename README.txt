#--------------------------------------------------------------#
#--------------------------------------------------------------#
# Pipeline for RNA-seq data analysis from raw fastq files      #
#--------------------------------------------------------------#
#--------------------------------------------------------------#

Author : Amandine Velt

==========================================================================
TABLE OF CONTENTS
==========================================================================

1. ORGANISATION OF THE PIPELINE SOURCES

2. REQUIREMENTS
  a. Software
  b. R and required libraries
  c. Data

3. INPUT

4. OUTPUT

5. USAGE

6. ORGANISATION OF THE SCRIPTS (which script runs what)

7. TO DO

==========================================================================

1. ORGANISATION OF THE PIPELINE SOURCES
=======================================
- sources

|---- README.txt
|---- config.sh
|---- run.sh
|---- src (dir)
|----|---- Alignment (dir)
|----|----|---- TophatInfos.R
|----|----|---- transcriptomeIndexVerification.sh
|----|---- Quality (dir)
|----|----|---- DrawQualityPlotsFastQC0.9.5.R
|----|----|---- GeneBodyCoverage.sh
|----|----|---- InnerDistance.sh
|----|----|---- OverallStatistics.sh
|----|----|---- QualityAnalysis.sh
|----|----|---- QualityInfos.R
|----|----|---- ReadsDistribution.sh
|----|----|---- StrandSpecificity.sh
|----|---- Quantification (dir)
|----|----|---- HTseqResultsAnalysis.sh
|----|---- RNAseqPipeline1Sample (dir)
|----|----|---- RNAseqPipeline1sample.sh
|----|----|---- RNAseqPipeline1sample_WithoutSpikes.sh
|----|---- RNAseqReport (dir)
|----|----|---- RNAseqReport.sh
|----|----|---- Files (dir)
|----|----|----|---- bibliography.bib
|----|----|----|---- design_file.txt
|----|----|----|---- genome_conversion.txt
|----|----|----|---- host_conversion.txt
|----|----|----|---- logo.jpg
|----|----|----|---- principal_file.tex
|----|---- Spikes (dir)
|----|----|---- FoldChange_graphic_script.R
|----|----|---- NRN_graphic_script.R
|----|----|---- SpikesAllSamplesAnalysis.sh
|----|----|---- SpikesBodyCoverage.sh
|----|----|---- Files (dir)
|----|----|----|---- Concentration.txt
|----|----|----|---- Expected_FoldChange.txt
|----|----|----|---- SAM_header_spikes.txt
|----|----|----|---- Spikes_Chromosomes_Length.txt
|----|----|----|---- Spikes_Length.txt
|----|----|----|---- Spikes_Length_sort.txt
|----|----|---- Report (dir)
|----|----|----|---- principal_file.tex
|----|----|----|---- SpikesReport.sh
|----|---- StatisticalAnalysis (dir)
|----|----|---- comp2cond.R
|----|----|---- comp2condDESeq2_withoutRep.R
|----|----|---- comp2cond_DESeq2.R
|----|----|---- comp2cond_var2samples.R
|----|----|---- RNAseqDataClustering.R
|----|----|---- RNAseqDataExploration.R
|----|----|---- RNAseqDataFormatingEnsVersion.R
|----|----|---- RNAseqDataStatistics.R
|----|----|---- RNAseqDataStatistics_DESeq2.R
|----|---- Utilities (dir)
|----|----|---- generateWig.sh
|----|----|---- gtf2bed.pl
|----|----|---- normalizeWig.pl
|----|----|---- sere.R
|----|----|---- utils.sh

2. REQUIREMENTS
===============

a. Software
-----------

- FastQC
Function : A quality control tool for high throughput sequence data
Download : http://www.bioinformatics.babraham.ac.uk/projects/download.html

- TopHat
Function : TopHat is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq readsto mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons. 
Download : TopHat or TopHat2 => http://ccb.jhu.edu/software/tophat/downloads/

- Bowtie
Function : Bowtie is an ultrafast, memory-efficient short read aligner. It aligns short DNA sequences (reads) to the human genome at a rate of over 25 million 35-bp reads per hour. Bowtie indexes the genome with a Burrows-Wheeler index to keep its memory footprint small.
Download : Bowtie => http://sourceforge.net/projects/bowtie-bio/files/bowtie/ or Bowtie 2 => http://sourceforge.net/projects/bowtie-bio/files/bowtie2/

- SAMtools
Function : SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format. 
Download : http://sourceforge.net/projects/samtools/files/samtools/

- HTSeq
Function : HTSeq is a Python package that provides infrastructure to process data from high-throughput sequencing assays. Given a SAM file with alignments and a GFF file with genomic features, this package counts how many reads map to each feature.
Download : https://pypi.python.org/pypi/HTSeq

- RSeQC
Function : RSeQC package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data.
Download : http://rseqc.sourceforge.net/#download-rseqc

- IGVTools
Function : The igvtools utility provides a set of tools for pre-processing data files. It can generated .wig and .tdf files.
Download : http://www.broadinstitute.org/software/igv/download (look at the "igvtools" part)

b. R and required libraries
-----------

- Biobase
Function : Base functions for Bioconductor
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("Biobase")

- DESeq
Function : Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("DESeq")

- DESeq2
Function : Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")

- edgeR
Function : Empirical analysis of digital gene expression data in R
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("edgeR")

- ade4
Function : Multivariate data analysis and graphical display
Installation with : install.packages("ade4")

- biomaRt
Function : Interface to BioMart databases (e.g. Ensembl, COSMIC ,Wormbase and Gramene)
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("biomaRt")

- GenomicRanges
Function : Representation and manipulation of genomic intervals
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges")

- rtracklayer
Function : R interface to genome browsers and their annotation tracks
Installation with : source("http://bioconductor.org/biocLite.R"); biocLite("rtracklayer")

- RColorBrewer
Function : Provides color schemes for maps (and other graphics)
Installation with : install.packages("RColorBrewer") 

- gplots
Function : Various R programming tools for plotting data
Installation with : install.packages('gplots')

c. Data
-------

- Genome 
|---- in fasta (in preference, downloaded from http://www.ensembl.org/info/data/ftp/index.html)
|---- bowtie or bowtie2 indexes (generated from the fasta downloaded before)
|---- gtf file (with genes/transcripts annotation. In preference, downloaded from http://www.ensembl.org/info/data/ftp/index.html)
- Spikes
|---- in fasta (downloaded from https://tools.lifetechnologies.com/content/sfs/manuals/ERCC92.zip)
|---- bowtie or bowtie2 indexes (generated from the fasta downloaded before)
|---- gtf file (https://tools.lifetechnologies.com/content/sfs/manuals/ERCC92.zip)

3. Input
========
Input data should be organized as such:

Path/to/Project_SNNNNN (outdir parameter of run.sh script)
|---- rawdata (dir)
|----|---- SampleID_1.fastq.gz
|----|---- SampleID_2.fastq.gz
|----|---- SampleID_3.fastq.gz
|----|---- SampleID_4.fastq.gz
|----|---- ...
|---- config.sh
|---- design_file.txt

Path/to/Project_SNNNNN is the only parameter of the script run.sh to run the RNAseqPipeline.
It must contains one rawdata folder with all the raw fastq files.
It must also contains the config file fills properly.
And finally, it must contains a design file of the experiment, like that :

#####-------------------------------------------------------
##### Sample_ID Sample_name Condition Contrast1  Contrast2
##### Sample1 Sample name1  Replicat1 1 1
##### Sample2 Sample name2  Replicat1 1 NA
##### Sample3 Sample name3  Replicat2 -1  -1
##### Sample4 Sample name4  Replicat2 -1  NA
#####-------------------------------------------------------

#####Contrast1 => comparison with replicates
#####Contrast2 => comparison without considering replicates
#####Separator : \t
#####If library type is paired-end, a line of the design file correspond to a pair of fastq (eg in design file, # # the sample_ID column of ${samplename}.R1.fastq.gz and ${samplename}.R2.fastq.gz is ${samplename}).
#####If library contains spikes, add "\_without_Spikes" at the end of sample_ID (indeed, spike's sequences are removed from the raw fastq file and the new fastq file, without spikes, is called ${sample_ID}_without_Spikes.fastq)

4. OUTPUT
=========

1. Output architecture
----------------------

Path/to/Project_SNNNNN/date (when you run RNAseqPipeline, a folder with the current date is created, this is the output folder of the run.)
|---- Alignment (dir)
|----|---- SampleID_1(_without_Spikes) (dir) (_without_Spikes is automatically added to the samplename if the raw fastq file contains spikes, because we remove spikes' sequences from the raw fastq)
|----|---- SampleID_2(_without_Spikes) (dir)
|----|---- ...
|---- Logs (dir)
|---- Processed_data (dir) (generated only if there are spikes in the experiment. The samples' raw sequences are seperated from spikes' sequences in this folder.)
|---- Quality (dir)
|----|---- fastqc (dir)
|----|---- RSeQC (dir)
|----|---- Spikes (dir)
|---- Quantification (dir)
|---- Report (dir)
|----|---- figures (dir)
|----|---- files (dir)
|----|---- Spikes (dir)
|---- Scripts (dir)
|----|---- src (dir)


2. Important output files
-------------------------

Path/to/Project_SNNNNN/date (when you run RNAseqPipeline, a folder with the current date is created, this is the output folder of the run.)
|---- Alignment (dir)
|----|---- SampleID_1(_without_Spikes) (dir) (_without_Spikes is automatically added to the samplename if the raw fastq file contains spikes, because we remove spikes' sequences from the raw fastq)
|----|----|---- SampleID_1(_without_Spikes)_alignment.bam
|----|----|---- SampleID_1(_without_Spikes)_alignment.bam.bai
|----|----|---- SampleID_1(_without_Spikes)_uniq_coverage_nonNormalized.wig.gz
|----|----|---- SampleID_1(_without_Spikes)_uniq_coverage_nonNormalized.tdf.gz
|----|---- SampleID_2(_without_Spikes) (dir)
|----|----|---- SampleID_2(_without_Spikes)_alignment.bam
|----|----|---- SampleID_2(_without_Spikes)_alignment.bam.bai
|----|----|---- SampleID_2(_without_Spikes)_uniq_coverage_nonNormalized.wig.gz
|----|----|---- SampleID_2(_without_Spikes)_uniq_coverage_nonNormalized.tdf.gz
|----|---- ...
|---- Logs (dir)
|---- Processed_data (dir) (generated only if there are spikes in the experiment. The samples' raw sequences are seperated from spikes' sequences in this folder.)
|---- Quality (dir)
|----|---- fastqc (dir)
|----|---- RSeQC (dir)
|----|---- Spikes (dir)
|----|---- all_statistics_results.txt
|---- Quantification (dir)
|----|---- htseq_combined.txt
|---- Report (dir)
|----|---- figures (dir)
|----|---- files (dir)
|----|---- SNNNNN_RNAseqDataExploration.pdf
|----|---- SNNNNN_report.pdf
|----|---- Spikes (dir)
|----|----|---- SNNNNN_spikes_report.pdf
|---- Scripts (dir)
|----|---- src (dir)

3. How the results/statistics have been done in all_statistics_results.txt
-------------------------------------------------------------------
/!\ this part of the README is valid with the use of TopHat2. The statistics of alignment are slightly different with TopHat1.

  a. for a single-read experiment

  - Total number of sequences: we count all the reads' header in the rawdata fastq files (e.g the number of "@..." in the fastq file)
  - With Spikes, we add the total number of sequences without spikes. We count all the reads' header in the processed fastq files (e.g the number of "@..." in the fastq file without spikes)
  - Unique number of sequences : we recover all the reads in the raw fastq file (or processed fastq file if there are spikes), we sort this file and we check the number of read's sequence present more
    than one time.
  - Total/Unique number of sequences : Total number of sequences divided by Unique number of sequences
  - With Spikes, we add the % of Spikes : we calculate the number of spikes (Total number of sequences - total number of sequences without spikes) and then we do -> number of Spikes * 100 / Total number of     sequences
  - Number of reads filtered out by Tophat : in the file "prep_reads.log" generated by TopHat, we recover the number of reads filtered out by the tool.
  - % of reads filtered out by Tophat : the percent is calculated compared to the total number of sequences or the total number of sequences without spikes, if we have spikes in the experiment.
  - Number of aligned reads : we recover the "Mapped" reads in the align_summary.txt file generated by TopHat2 and we calculate the % compared to the total number of sequences or the total number of sequences without spikes, if we have spikes in the experiment.
  - Number of multiple aligned reads : we recover the multiple aligned reads in the align_summary.txt file generated by TopHat2 and we calculate the % compared to the total number of aligned reads.
  - Number of uniquely aligned reads : we calculate the number of uniquely aligned reads (Number of aligned reads - Number of multiple aligned reads) and we calculate the % compared to the total number of aligned reads.
  - The density (Tags/Kb) in exons, introns and TSS up 10kb / TES down 10kb are calculated by the tool RSeQC, on the reads of spikes if there are spikes, else on the total reads.
  - Number of assigned reads, Number of no feature reads, Number of ambiguous reads and Number of multiple alignments are calculated by the tool HTSeq and we calculate the % compared to the number of uniquely aligned reads (because HTSeq works on these reads only)

  b. for a paired-end experiment, there are some supplementary columns

  - Number of aligned reads calculated by TopHat2 : we recover the "Aligned pairs" statistic in the align_summary.txt file generated by TopHat2 and we calculate the % compared to the total number of sequences or the total number of sequences without spikes, if we have spikes in the experiment.
  - Number of multiple aligned reads calculated by TopHat2 : we recover the multiple aligned reads in the align_summary.txt file from the "aligned pairs" generated by TopHat2 and we calculate the % compared to the total number of aligned reads calculated by TopHat2.
  - Number of uniquely aligned reads : we calculate the number of uniquely aligned reads (Number of aligned reads calculated by TopHat2 - Number of multiple aligned reads calculated by TopHat2) and we calculate the % compared to the total number of aligned reads calculated by TopHat2.
  - Number of discordant pairs alignment calculated by TopHat2 : we recover the "discordant" statistic in the align_summary.txt file generated by TopHat2 and we calculate the % compared to the total number of aligned reads calculated by TopHat2.
  

  - Number of aligned reads : because TopHat2 did his statistics on the "aligned pairs", it underestimates the number of aligned reads. So we count all the reads which are mapped on the genome (with its pair or not) and we calculate the % compared to the total number of sequences or the total number of sequences without spikes, if we have spikes in the experiment.
  - Number of uniquely aligned reads : because TopHat2 did his statistics on the "aligned pairs", it underestimates the number of uniquely aligned reads. So we count all the reads which are mapped only one time on the genome (with its pair or not) and we calculate the % compared to the total number of aligned reads.
  - Number of multiple aligned reads : we calculate the number of of multiple aligned reads (Number of aligned reads - Number of uniquely aligned reads) and we calculate the % compared to the total number of aligned reads.

5. USAGE
========

/ !!! \ run the RNAseqPipeline with slurm.

outdir_run="PATH/TO/outdir"
cd $outdir_run
source config.sh

'sbatch --job-name "RNAseqPipeline" -p $SLURM_PARTITION -o "$outdir_run/nohup.allsamples" <<-EOF
#!/bin/bash
bash $SCRIPTSDIR/run.sh -o "$outdir_run"
EOF'

- outdir is the full path to the folder containing the three essential folders/files:
|----the rawdata folder containing all the raw fastq files
|----the config file fills properly
|----the design file of the experiment, fills properly

6. ORGANISATION OF THE SCRIPTS (which script runs what)
=======================================================
|-- run.sh
  |---- src/RNAseqPipeline1Sample/RNAseqPipeline1sample.sh
  |----|---- src/Quality/ReadsDistribution.sh
  |----|---- src/Quality/StrandSpecificity.sh
  |----|---- src/Quality/GeneBodyCoverage.sh
  |----|---- src/Quality/InnerDistance.sh (if paired-end data)
  |----|---- src/Quality/OverallStatistics.sh
  |---- src/Quality/QualityAnalysis.sh
  |----|---- src/Quality/DrawQualityPlotsFastQC0.9.5.R
  |---- src/Spikes/SpikesAllSamplesAnalysis.sh
  |----|---- src/Spikes/SpikesBodyCoverage.sh
  |----|---- src/Spikes/NRN_graphic_script.R
  |----|---- src/Spikes/FoldChange_graphic_script.R
  |----|---- src/Spikes/Report/SpikesReport.sh
  |---- src/Quality/QualityInfos.R
  |---- src/Alignment/TophatInfos.R
  |---- src/Quantification/HTseqResultsAnalysis.sh
  |---- src/RNAseqReport/RNAseqReport.sh
  |----|---- src/StatisticalAnalysis/RNAseqDataClustering.R
  |----|---- src/StatisticalAnalysis/RNAseqDataFormatingEnsVersion.R
  |----|---- src/StatisticalAnalysis/RNAseqDataExploration.R
  |----|---- src/Utilities/sere.R
  |----|---- src/StatisticalAnalysis/RNAseqDataStatistics.R
  |----|----|---- src/StatisticalAnalysis/comp2cond.R OR src/StatisticalAnalysis/comp2cond_var2samples.R
  |----|---- src/StatisticalAnalysis/RNAseqDataStatistics_DESeq2.R
  |----|----|---- src/StatisticalAnalysis/comp2cond_DESeq2.R OR src/StatisticalAnalysis/comp2condDESeq2_withoutRep.R
  |---- src/Utilities/utils.sh && wigTdfGeneration
  |----|---- src/Utilities/generateWig.sh
  |----|----|---- src/Utilities/normalizeWig.pl

7. TO DO
========

- Make transcripts analysis