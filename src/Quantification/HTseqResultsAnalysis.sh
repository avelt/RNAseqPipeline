#! /bin/sh

#Author : keime

###############################################################################################
#This script prints several statistics on HTseq results in $outdir/htseq_analysis/htseq_stats_summary.txt
#and creates a file with all HTSeq results combined in $outdir/htseq_analysis/htseq_combined.txt
###############################################################################################

outdir=$1

if [ "$outdir" == "" ]
then
        echo "Usage : sh HTseqResultsAnalysis.sh path/to/outdir"
        exit    
fi

#outdir : path/to/directory where there is a Alignment directory 
#(in the Alignment directory there should be a directory for each sample, with tophat results in this directory)


#create output directory
cd $outdir

#ceate a file to output several statistics
echo "Sample name	Number of assigned reads	Number of no feature reads	Number of ambiguous reads	Number of multiple alignments" > $outdir/Quantification/htseq_stats_summary.txt

#create a title file
printf "Ensembl gene id" > tmptitre

#for each directory in the Alignment directory 
for dir in `ls -d $outdir/Alignment/*/`
do
	echo "Analysis of directory $dir"
	#name of the directory
	name=`echo $dir | perl -e '$in=<>;@s=split("/",$in);print "$s[$#split-1]";'`

	if [ -f $dir/"$name"_htseq.txt ]
	then

	printf "$name" >> $outdir/Quantification/htseq_stats_summary.txt
	
	#print the name into a file to create the title of the combined htseq file
	printf "	$name" >> tmptitre 
	
	#create a file with gene names
	grep 'ENS' $dir/"$name"_htseq.txt | cut -f 1 > tmpgenes_"$name"
	
	#create a file with corresponding occurrence number
	grep 'ENS' $dir/"$name"_htseq.txt | cut -f 2  > tmpocc_"$name"
	
	
	assigned=`grep 'ENS' $dir/"$name"_htseq.txt | cut -f 2 | perl -e '$sum=0;while ($line=<>) {$sum=$sum+$line;}; print $sum; '`
	nofeature=`grep no_feature $dir/"$name"_htseq.txt | cut -f 2`
	ambiguous=`grep ambiguous $dir/"$name"_htseq.txt | cut -f 2`
	notunique=`grep alignment_not_unique $dir/"$name"_htseq.txt | cut -f 2`
	echo "	$assigned	$nofeature	$ambiguous	$notunique" >> $outdir/Quantification/htseq_stats_summary.txt
	
	else

		echo "/!\ /!\ /!\ $dir/${name}_htseq.txt doesnt exist !! this sample is not added to the htseq_combined.txt file ! /!\ /!\ /!\ "

	fi
done

#compare the name of all genes to verify that there are the same
printf "Check if all files contain the same genes identification numbers..."
comp=`ls tmpgenes* | head -1` #le premier fichier, auquel seront compares tous les autres
for i in `ls tmpgenes*`
do 
	diff $i $comp >> test
done


#if all genes are the same, combine all HTseq results file 
if (file test | grep -v empty) >> /dev/null
	then
		echo "The identification numbers are not the same in all files !!!"
	else
		echo " OK"
		echo "Merge all files"
		paste $comp tmpocc_* > tmpcombined
		echo >> tmptitre
		cat tmptitre tmpcombined > $outdir/Quantification/htseq_combined.txt
fi

#remove temporary files
rm -f tmp* test

