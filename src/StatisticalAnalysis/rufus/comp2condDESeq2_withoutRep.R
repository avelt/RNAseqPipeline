#Author : Keime, Velt
#Date : 12/2013

#function to compare 2 conditions
#comparisons are done only on samples compared, and a new normalization is done
#and does not match that normalization values of alldata.

comp2cond_DESeq2_withoutRep = function(ddsComp, c, tlfc, tpadj, genes.res, tlfc, tpadj){

	#name of the comparison (keep only condition1_vs_condition2, ie without the word Condition)
	cname = substr(c,11,nchar(c)) 
	#name of the compared conditions
	cond1 = strsplit(cname,"_vs_")[[1]][2]
	cond2 = strsplit(cname,"_vs_")[[1]][1]
					
	#results of this comparison
	res=results(ddsComp, c)	
					
	#MAplot
	png(paste("figures/MA_",c,"_DESeq2.png",sep=""))
	#changing default MAplot to add a threshold on log2 fold-change
	col=ifelse(mcols(ddsComp)[, paste("WaldAdjPvalue_",c,sep="")] < tpadj & abs(mcols(ddsComp)[, c]) > tlfc, "red", "black")
	suppressWarnings(plotMA(ddsComp, c, col=col))
	dev.off()

	#p-value histogram
	png(paste("figures/histpval_",c,"_DESeq2.png",sep="")) 
	suppressWarnings(hist(res$pvalue,100, main="",xlab="p-value"))
	dev.off()
						
	#print number of significantly differentially expressed genes
	cat(paste("\n\n \\textbf{Comparison between ",as.character(cond2), " and ", as.character(cond1), " with DESeq2.} \\newline", sep=""))
	cat(paste("\n Total number of significantly differentially expressed genes: ", sum(res$padj<tpadj & abs(res$log2FoldChange)>tlfc,na.rm=T), "\\newline"))
	cat(paste("\n Number of overexpressed genes: ", sum(res$padj<tpadj & res$log2FoldChange>tlfc,na.rm=T), "\\newline"))
	cat(paste("\n Number of underexpressed genes: ", sum(res$padj<tpadj & res$log2FoldChange<(-tlfc),na.rm=T), "\\newline"))
							
	#vector containing a summary of the results (-1 if the gene is down and with 1 if the gene is up)
	brief = data.frame(rep(0,length=dim(ddsComp)[1]), row.names=rownames(ddsComp))
	colnames(brief)[1] = cname
	brief[ rownames(brief) %in% rownames(res)[res$padj<tpadj & !is.na(res$padj) & res$log2FoldChange>tlfc] ,1 ] = 1
	brief[ rownames(brief) %in% rownames(res)[res$padj<tpadj & !is.na(res$padj) & res$log2FoldChange<(-tlfc)], 1 ] = -1
							
	#append the results to the alldata data.frame
	all = merge(genes.res, as.data.frame(res[,c(-1)]), by.x="Ensembl gene id", by.y="row.names", all=T)
	all = all[order(all$padj),]
						
	nbcol = dim(all)[2]
					
	colnames(all)[(nbcol-3):nbcol] = c(paste("log2(",cond2,"/",cond1,")_withoutRep", sep=""), paste("standard error for the log2(",cond2,"/",cond1,")_withoutRep", sep=""), paste("P-value (",cond2," vs ",cond1,")_withoutRep", sep=""), paste("Adjusted p-value ", "(",cond2," vs ",cond1,")_withoutRep", sep=""))
	return(all)
	
}
