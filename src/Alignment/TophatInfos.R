#Author : keime

###############################################################################################
#This function creates a tab delimited tophat_summary.file.txt in dir directory
#with information on tophat alignment :
#total nb of reads
#nb of reads filtered out by tophat
#% of reads filtered out by tophat
#nb of aligned reads
#% of aligned reads
#nb of uniquely aligned reads
#% of uniquely aligned reads
#nb of multiple aligned reads
#% of multiple aligned reads

#This script uses information from dir/tophat_stats.txt file 
#################################################################################################



tophatinfos = function(dir,paired,tophat2){
	#read stats.txt file from dir directory
	data = read.table(paste(dir,"/tophat_stats.txt",sep=""), sep="|")
	
	nblines = dim(data)[1] #nb of lines of the data data.frame

	# if paired-end data, we have a supplementary information, which is : discordant pairs alignment
	if ((paired=="yes") && (tophat2=="yes")){

		nbinfos = 9 #nb of information per sample
		#infos data.frame with information on the alignment
		infos = data.frame(
		data$V2[seq(1,nblines,nbinfos)], #directory => same that below
      	data$V2[seq(2,nblines,nbinfos)], #nb of reads => same that below
		as.numeric(as.character(data$V2[seq(3,nblines,nbinfos)])), #nb of aligned reads calculated by tophat  => add by use of tophat2 only
		as.numeric(as.character(data$V2[seq(4,nblines,nbinfos)])), #nb of uniquely aligned reads calculated by tophat  => add by use of tophat2 only
		as.numeric(as.character(data$V2[seq(5,nblines,nbinfos)])), #nb of multiple aligned reads calculated by tophat  => add by use of tophat2 only
		as.numeric(as.character(data$V2[seq(6,nblines,nbinfos)])), #nb of aligned reads => same that below
		as.numeric(as.character(data$V2[seq(7,nblines,nbinfos)])), #nb of uniquely aligned reads => same that below
		as.numeric(as.character(data$V2[seq(8,nblines,nbinfos)])), #nb of multiple aligned reads => same that below
		as.numeric(as.character(data$V2[seq(9,nblines,nbinfos)])) #nb of discordant pairs => add by use of tophat2 only
	)

	} else {

		nbinfos = 5 #nb of information per sample
		infos = data.frame(
		data$V2[seq(1,nblines,nbinfos)], #directory
      	data$V2[seq(2,nblines,nbinfos)], #nb of reads
		as.numeric(as.character(data$V2[seq(3,nblines,nbinfos)])), #nb of aligned reads
		as.numeric(as.character(data$V2[seq(4,nblines,nbinfos)])), #nb of uniquely aligned reads
		as.numeric(as.character(data$V2[seq(5,nblines,nbinfos)])) #nb of multiple aligned reads
		)
	}
	
	#keep only sample ID (not entire path)
	pathsplited=strsplit(as.character(infos[,1]),"/")
	lastelt=function(x){x[length(x)]}
	samplesid=sapply(pathsplited, lastelt)
	
	#nb reads filtered out by tophat
	readssplited=strsplit(as.character(infos[,2]),"out of")
	firsteltnum=function(x){as.numeric(x[1])}
	filteredreads=sapply(readssplited, firsteltnum)
    
	#total nb of reads
	tmptotalreads=sapply(readssplited, lastelt)
	tmptotalreadssplited = strsplit(tmptotalreads, "reads")
	totalreads = sapply(tmptotalreadssplited, firsteltnum)
  
	#calculate % of filtered out, aligned, uniquely aligned and multiple aligned reads
	pfilter = round(filteredreads*100/totalreads, digits=2)
	palign = round(infos[,3]*100/totalreads, digits=2)
	pualign = round(infos[,4]*100/totalreads, digits=2)
	pmalign = round(infos[,5]*100/totalreads, digits=2) 
	if ((paired=="yes") && (tophat2=="yes")){
		pdpairs = round(infos[,9]*100/totalreads, digits=2)
		palign_2 = round(infos[,6]*100/totalreads, digits=2)
		pualign_2 = round(infos[,7]*100/totalreads, digits=2)
		pmalign_2 = round(infos[,8]*100/totalreads, digits=2)
	}         
  
  #all information in a data.frame
  	if ((paired=="yes") && (tophat2=="yes")){

		infos.ok = data.frame(samplesid, totalreads, filteredreads, pfilter, infos[,3], palign, infos[,4], pualign, infos[,5], pmalign, infos[,9], pdpairs, infos[,6], palign_2, infos[,7], pualign_2, infos[,8], pmalign_2)
		colnames(infos.ok) = c("Sample ID", "Total number of reads", "Number of reads filtered out by Tophat", "% of reads filtered out by Tophat",
			"Number of aligned reads calculated by TopHat2", "% of aligned reads calculated by TopHat2", 
			"Number of uniquely aligned reads calculated by TopHat2", "% of uniquely aligned reads calculated by TopHat2",
			"Number of multiple aligned reads calculated by TopHat2", "% of multiple aligned reads calculated by TopHat2",
			"Number of discordant pairs alignment calculated by TopHat2", "% of discordant pairs alignment calculated by TopHat2",
			"Number of aligned reads", "% of aligned reads", 
			"Number of uniquely aligned reads", "% of uniquely aligned reads",
			"Number of multiple aligned reads", "% of multiple aligned reads")   

	} else {

		infos.ok = data.frame(samplesid, totalreads, filteredreads, pfilter, infos[,3], palign, infos[,4], pualign, infos[,5], pmalign)
		colnames(infos.ok) = c("Sample ID", "Total number of reads", "Number of reads filtered out by Tophat", "% of reads filtered out by Tophat",
			"Number of aligned reads", "% of aligned reads", 
			"Number of uniquely aligned reads", "% of uniquely aligned reads",
			"Number of multiple aligned reads", "% of multiple aligned reads"
		)                     

	}
  
	
	#write infos data.frame in a file dir/stats_summary.txt
	write.table(infos.ok, file=paste(dir,"/tophat_stats_summary.txt",sep=""), quote=F, sep="\t", row.names=F, dec=",")
}



