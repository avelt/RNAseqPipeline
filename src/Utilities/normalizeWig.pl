#! /usr/bin/perl

#Author : keime
#Date : 2013/09

#########################################################################################
#This script normalize wig files by dividing each value by the normalization factor provided
#and exclude spikes sequences from wig file (beginning by ERCC)
#Usage : perl normalizeWig file.wig normfactor trackline > file_norm.wig
#file.wig : a wig file to be normalized
#normfactor : the normalization factor, eg. provided by DESeq sizeFactors() function
#trackline : characters that will be printed in the track definition line
#file_norm.wig output file (same format as input, all values in column 4 are divided by the size factor
#########################################################################################

#open wig file
open(WIG,"$ARGV[0]") or die ("File $ARGV[0] can not be opened");

#normalization factor
$normfac=$ARGV[1];

#read first line of wig file
$line = <WIG>; 

#print track definition line provided in input
print $ARGV[2]."\n";

while( defined( $line = <WIG> ) ) { #for each line of wig file
	if ($line !~ /^ERCC/){ #if it does not correspond to a spike sequence
	
		if ($line eq /variableStep.*/){ #the variableStep lines are kept because declaration lines
		
			print "$line\n";
		
		} else {
	
			@linesplit = split("\t",$line);
			print "$linesplit[0]\t"; #print first columns as in input file
			$normcov = $linesplit[1]/$normfac; #normalized coverage
			print "$normcov\n"; #2nd column : input value / normalization factor
		
		}
	}
}

#close wig file
close(WIG);
