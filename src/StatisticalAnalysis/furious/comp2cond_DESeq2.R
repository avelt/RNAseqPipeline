#Author : Keime, Velt
#Date : 12/2013

#function to compare 2 conditions
#comparisons are done only on samples compared, and a new normalization is done
#and does not match that normalization values of alldata.

comp2cond_DESeq2 = function(rawdataOK, dds, Cond1, Cond2, tlfc, tpadj, genes.res, counter, index_cond1_for_scatterplot, index_cond2_for_scatterplot){
	
	if (counter > 1){

	#results of this comparison
	compName=paste(Cond2,paste("_vs_",Cond1, "_no_replicates",sep=""),sep="")
	res=results(dds, contrast=c("Condition", Cond2,Cond1))

	#MAplot
	png(paste("figures/MA_",compName,"_DESeq2.png",sep=""))
	plotMA(res, alpha=tpadj, main="DESeq2", ylim=c(-2,2))
	dev.off()

	#p-value histogram
	png(paste("figures/histpval_",compName,"_DESeq2.png",sep="")) 
	suppressWarnings(hist(res$pvalue,100, main="",xlab="p-value"))
	dev.off()

	#plot of counts for gene
	png(paste("figures/plot_of_counts_gene_for_",compName,"_DESeq2.png",sep=""))
	plotCounts(dds,gene=which.min(res$padj), intgroup="Condition")
	dev.off()

	#print number of significantly differentially expressed genes
	cat(paste("\n\n \\textbf{Comparison between ",as.character(Cond2), " and ", as.character(Cond1), " with DESeq2.} \\newline", sep=""))
	cat(paste("\n Total number of significantly differentially expressed genes: ", sum(res$padj<tpadj & abs(res$log2FoldChange)>tlfc,na.rm=T), "\\newline"))
	cat(paste("\n Number of overexpressed genes: ", sum(res$padj<tpadj & res$log2FoldChange>tlfc,na.rm=T), "\\newline"))
	cat(paste("\n Number of underexpressed genes: ", sum(res$padj<tpadj & res$log2FoldChange<(-tlfc),na.rm=T), "\\newline"))

	countCond1=rowMeans(rawdataOK[,c(index_cond1_for_scatterplot)])
	countCond2=rowMeans(rawdataOK[,c(index_cond2_for_scatterplot)])

	#Recovery of the greatest value for scatterplot
	MaxValue1=max(countCond1)
	#we select res$baseMean > 1 because if res$baseMean<1, log(res$baseMean) <0 et xlim-ylim <0
	MinValue1=min(countCond1[countCond1>=1])
	MaxValue2=max(countCond2)
	MinValue2=min(countCond2[countCond2>=1])
	
	if (MinValue1 <= MinValue2) {MinValue=MinValue1 } else { MinValue=MinValue2}
	if (MaxValue1 >= MaxValue2) {MaxValue=MaxValue1 } else { MaxValue=MaxValue2}

	png(paste("figures/",Cond2,"_vs_",Cond1,"_scatterplot.png", sep=""))
	suppressWarnings(plot(countCond1, countCond2, log="xy", pch=20, cex=0.6, col=ifelse(res$padj<tpadj & abs(res$log2FoldChange) > tlfc, "red", "black") , cex.lab=1, cex.axis=1, xlab=Cond1, ylab=Cond2, xlim = c(MinValue,MaxValue), ylim = c(MinValue,MaxValue)))
	legend("bottomright", col=c(2,1), legend=c("FDR<threshold  and |lfc|>threshold ","FDR>threshold  or |lfc|<threshold "), pch=20, cex=0.9)
	abline(0, 1)
	dev.off()

	} else {

	#results of this comparison
	compName=paste(Cond2,paste("_vs_",Cond1,sep=""),sep="")
	res=results(dds, contrast=c("Condition", Cond2,Cond1))

	#MAplot
	png(paste("figures/MA_",compName,"_no_replicates_DESeq2.png",sep=""))
	plotMA(res, alpha=tpadj, main="DESeq2", ylim=c(-2,2))
	dev.off()

	#p-value histogram
	png(paste("figures/histpval_",compName,"_no_replicates_DESeq2.png",sep="")) 
	suppressWarnings(hist(res$pvalue,100, main="",xlab="p-value"))
	dev.off()

	#plot of counts for gene
	png(paste("figures/plot_of_counts_gene_for_",compName,"_no_replicates_DESeq2.png",sep=""))
	plotCounts(dds,gene=which.min(res$padj), intgroup="Condition")
	dev.off()

	#print number of significantly differentially expressed genes
	cat(paste("\n\n \\textbf{Comparison between ",as.character(Cond2), " and ", as.character(Cond1), " with DESeq2. !! WARNING : no replicates used !!} \\newline", sep=""))
	cat(paste("\n Total number of significantly differentially expressed genes: ", sum(res$padj<tpadj & abs(res$log2FoldChange)>tlfc,na.rm=T), "\\newline"))
	cat(paste("\n Number of overexpressed genes: ", sum(res$padj<tpadj & res$log2FoldChange>tlfc,na.rm=T), "\\newline"))
	cat(paste("\n Number of underexpressed genes: ", sum(res$padj<tpadj & res$log2FoldChange<(-tlfc),na.rm=T), "\\newline"))

	countCond1=rawdataOK[,c(index_cond1_for_scatterplot)]
	countCond2=rawdataOK[,c(index_cond2_for_scatterplot)]

	#Recovery of the greatest value for scatterplot
	MaxValue1=max(countCond1)
	#we select res$baseMean > 1 because if res$baseMean<1, log(res$baseMean) <0 et xlim-ylim <0
	MinValue1=min(countCond1[countCond1>=1])
	MaxValue2=max(countCond2)
	MinValue2=min(countCond2[countCond2>=1])
	
	if (MinValue1 <= MinValue2) {MinValue=MinValue1 } else { MinValue=MinValue2}
	if (MaxValue1 >= MaxValue2) {MaxValue=MaxValue1 } else { MaxValue=MaxValue2}

	png(paste("figures/",Cond2,"_vs_",Cond1,"_scatterplot_no_replicates.png", sep=""))
	suppressWarnings(plot(countCond1, countCond2, log="xy", pch=20, cex=0.6, col=ifelse(res$padj<tpadj & abs(res$log2FoldChange) > tlfc, "red", "black") , cex.lab=1, cex.axis=1, xlab=Cond1, ylab=Cond2, xlim = c(MinValue,MaxValue), ylim = c(MinValue,MaxValue)))
	legend("bottomright", col=c(2,1), legend=c("FDR<threshold and |lfc|>threshold ","FDR>threshold  or |lfc|<threshold"), pch=20, cex=0.9)
	abline(0, 1)
	dev.off()

	}
						
	#vector containing a summary of the results (-1 if the gene is down and with 1 if the gene is up)
	brief = data.frame(rep(0,length=dim(dds)[1]), row.names=rownames(dds))
	colnames(brief)[1] = compName
	brief[ rownames(brief) %in% rownames(res)[res$padj<tpadj & !is.na(res$padj) & res$log2FoldChange>tlfc] ,1 ] = 1
	brief[ rownames(brief) %in% rownames(res)[res$padj<tpadj & !is.na(res$padj) & res$log2FoldChange<(-tlfc)], 1 ] = -1
							
	#append the results to the alldata data.frame
	all = merge(genes.res, as.data.frame(res[,c(-1,-4)]), by.x="Ensembl gene id", by.y="row.names", all=T)
	all = all[order(all$padj),]
						
	nbcol = dim(all)[2]
						
	colnames(all)[(nbcol-3):nbcol] = c(paste("log2(",Cond2,"/",Cond1,")", sep=""), paste("standard error for the log2(",Cond2,"/",Cond1,")", sep=""), paste("P-value (",Cond2," vs ",Cond1,")", sep=""), paste("Adjusted p-value ", "(",Cond2," vs ",Cond1,")", sep=""))
	return(list(allres=brief, genes.res=all))
	
}
