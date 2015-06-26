#Author : keime

###############################################################################################
#This function creates a tab delimited stats_summary.file.txt in dir directory
#with information on the samples and quality-related parameters :
#sequencer
#flow cell 
#sequence length
#line
#sample id
#total nb of reads
#unique nb of reads
#total /unique nb of reads

#This script uses information from dir/stats.txt file where the file name are *.fastq.gz
#################################################################################################



qualityinfos = function(dir, spikes){
	#read stats.txt file from dir directory
	data = read.table(paste(dir,"/stats.txt",sep=""), sep="|")

	# if spikes=yes, we have a supplementary information, which is : total reads without spikes
	if (spikes=="yes") {

		nblines = dim(data)[1] #nb of lines of the data data.frame
		nbinfos = 8 #nb of information per sample
		#infos data.frame with information on the samples and quality-related parameters
		infos = data.frame(
			data$V2[seq(1,nblines,nbinfos)], #sample ID
			data$V2[seq(2,nblines,nbinfos)], #sequencer
			data$V2[seq(3,nblines,nbinfos)], #flow cell
			data$V2[seq(4,nblines,nbinfos)], #line
			data$V2[seq(5,nblines,nbinfos)], #read length
			as.numeric(as.character(data$V2[seq(6,nblines,nbinfos)])), #total nb of reads
			as.numeric(as.character(data$V2[seq(7,nblines,nbinfos)])), #total nb of reads without spikes
			as.numeric(as.character(data$V2[seq(8,nblines,nbinfos)])) #unique nb of reads
		)

	 } else {

		nblines = dim(data)[1] #nb of lines of the data data.frame
		nbinfos = 7 #nb of information per sample
		#infos data.frame with information on the samples and quality-related parameters
		infos = data.frame(
			data$V2[seq(1,nblines,nbinfos)], #sample ID
			data$V2[seq(2,nblines,nbinfos)], #sequencer
			data$V2[seq(3,nblines,nbinfos)], #flow cell
			data$V2[seq(4,nblines,nbinfos)], #line
			data$V2[seq(5,nblines,nbinfos)], #read length
			as.numeric(as.character(data$V2[seq(6,nblines,nbinfos)])), #total nb of reads
			as.numeric(as.character(data$V2[seq(7,nblines,nbinfos)])) #unique nb of reads
		)
	}
	
	
	#keep only sample ID (not entire path)
	pathsplited=strsplit(as.character(infos[,1]),"/")
	lastelt=function(x){x[length(x)]}
	files=lapply(pathsplited, lastelt)
	beforefastq = function(x){unlist(strsplit(x, ".fastq.gz"))}
	infos[,1] = unlist(lapply(files, beforefastq))
	
	#calculate percent total/quality filtered sequences (ptq)
	#ptq = round(infos[,7]/infos[,6]*100, digits=2)
	
	#calculate total/unique nb of sequences
	if (spikes=="yes") {

		infos = data.frame(infos, round(infos[,7]/infos[,8], digits=2))
		colnames(infos) = c(as.character(data$V1[1:8]),"Total/Unique number of sequences")

	} else {
		infos = data.frame(infos, round(infos[,6]/infos[,7], digits=2))
		colnames(infos) = c(as.character(data$V1[1:7]),"Total/Unique number of sequences")
	}
	
	#write infos data.frame in a file dir/stats_summary.txt
	write.table(infos, file=paste(dir,"/stats_summary.txt",sep=""), quote=F, sep="\t", row.names=F)
}




