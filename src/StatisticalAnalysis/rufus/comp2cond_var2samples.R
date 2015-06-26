#Author : Keime, Velt
#Date : 02/2013

#function to compare 2 conditions
#estimate variance on the 2 samples considered only

#cds : a CountDataSet object
#i, j : sample numbers to compare (i correspond to condition cond1 and j to condition cond2)
#cond1 : name of the first condition to compare (this must match one of the conditions from phenoData(cds)@data$condition )
#cond2 : name of a second condition to compare with cond1 (this must match one of the conditions from phenoData(cds)@data$condition )
#allres : a data.frame with as many lines as the nb of genes (dim(cds)[1]) and as many columns as the nb of comparisons you want to perform
#allres4evenn : a data.frame with as many lines as the nb of genes (dim(cds)[1]) and nb of columns=(nb of comparisons you want to perform)x2+1
#nb : column nb in allres to store the corresponding results (1 if overexpressed, -1 if underexpressed)
#genes.res : data.frame on which to merge results of this analysis


comp2cond.var2samples = function(cds, i, j, cond1, cond2, allres, nb, genes.res) {
	
	cds = estimateDispersions( cds[,c(i,j)], method="blind", sharingMode="fit-only") #as we do not have replicates
	
	#Recovery of the greatest value of cds
		
	res=nbinomTest( cds, cond1, cond2 )
	res = res[! is.na(res$pval),]
	
	#Recovery of the greatest value for scatterplot
	MaxValue1=max(res$baseMeanA)
	#we select res$baseMean > 1 because if res$baseMean<1, log(res$baseMean) <0 et xlim-ylim <0
	MinValue1=min(res$baseMeanA[res$baseMeanA>=1])
	MaxValue2=max(res$baseMeanB)
	MinValue2=min(res$baseMeanB[res$baseMeanB>=1])
	
	if (MinValue1 <= MinValue2) {MinValue=MinValue1 } else { MinValue=MinValue2}
	if (MaxValue1 >= MaxValue2) {MaxValue=MaxValue1 } else { MaxValue=MaxValue2}
	
	#MinValueOk=log(MinValue)
	#MaxValueOk=log(MaxValue)
	
	#print nb of differentially expressed genes	
	cat(paste("\n\n \\textbf{Comparison between ",as.character(cond2), " and ", as.character(cond1), " with DESeq } \\newline", sep=""))
	cat(paste("\n Total number of significantly differentially expressed genes: ", sum(res$padj<0.05 & abs(res$log2FoldChange)>1  ), "\\newline" ))
	cat(paste("\n Number of overexpressed genes: ", sum(res$padj<0.05 & res$log2FoldChange>1 ), "\\newline "))
	cat(paste("\n Number of underexpressed genes: ", sum(res$padj<0.05 & res$log2FoldChange<(-1)), "\\newline "))
	
	
	#fill allres with the corresponding results
	allres[   which (rownames(allres) %in% res$id[res$padj<0.05 & res$log2FoldChange>1]) , nb ] = 1
	allres[   which (rownames(allres) %in% res$id[res$padj<0.05 & res$log2FoldChange<(-1)]) , nb ] = -1
	colnames(allres)[nb] = paste(cond2,"vs",cond1)
	
	#scatter plot
	pdf(paste("figures/",cond2,"_vs_",cond1,"_scatterplot.pdf", sep=""))
	suppressWarnings(plot(res$baseMeanA, res$baseMeanB, log="xy", pch=10, cex=0.6, col=ifelse(res$padj<0.05 & abs(res$log2FoldChange) > 1, "red", "black") , cex.lab=1, cex.axis=1, 
			xlab= cond1, ylab=cond2, xlim = c(MinValue,MaxValue), ylim = c(MinValue,MaxValue)))
	legend("bottomright", col=c(2,1), legend=c("FDR<0.05 and |lfc|>1","FDR>0.05 or |lfc|<1"), pch=20, cex=0.9)
	dev.off()
	
	png(paste("figures/",cond2,"_vs_",cond1,"_scatterplot.png", sep=""))
	suppressWarnings(plot(res$baseMeanA, res$baseMeanB, log="xy", pch=20, cex=0.6, col=ifelse(res$padj<0.05 & abs(res$log2FoldChange) > 1, "red", "black") , cex.lab=1, cex.axis=1, 
			xlab=cond1, ylab=cond2, xlim = c(MinValue,MaxValue), ylim = c(MinValue,MaxValue)))
	legend("bottomright", col=c(2,1), legend=c("FDR<0.05 and |lfc|>1","FDR>0.05 or |lfc|<1"), pch=20, cex=0.9)
	dev.off()
	
	#MA plot
	png(paste("figures/",cond2,"_vs_",cond1,"_MA.png", sep=""))
	suppressWarnings(plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=0.6, col=ifelse(res$padj<0.05 & abs(res$log2FoldChange) > 1, "red", "black") , 
			cex.lab=1, cex.axis=1, xlab="Mean of normalized counts", ylab="log2 fold change", ylim=c(-10,10)))
	legend("bottomright", col=c(2,1), legend=c("FDR<0.05 and |lfc|>1","FDR>0.05 or |lfc|<1"), pch=20, cex=0.9)
	abline(h=0)
	#abline(h=1, lty=2)
	#abline(h=-1, lty=2)
	dev.off()
		
	#histogram of p-values
	pdf(paste("figures/",cond2,"_vs_",cond1,"_histpval.pdf", sep=""))
	hist(res$pval)
	dev.off()
	
	#add the corresponding results to a data.frame
	colnames(res)[1] = "Ensembl gene id"
	genes.res.new = merge(genes.res, res, by.x="Ensembl gene id", all=T)
	nbcol = dim(genes.res.new)[2]
	genes.res.ok = genes.res.new[order(genes.res.new$padj),-c(nbcol-4, nbcol-5, nbcol-6)]
	nbcol = dim(genes.res.ok)[2]
	colnames(genes.res.ok)[(nbcol-3):nbcol] = c(paste(cond2,"/",cond1, sep=""), paste("log2(",cond2,"/",cond1,")", sep=""), paste("P-value (",cond2," vs ",cond1,")", sep=""), paste("Adjusted p-value ", "(",cond2," vs ",cond1,")", sep=""))
	genes.res = genes.res.ok
	
	return(list(allres=allres, genes.res=genes.res))
}
