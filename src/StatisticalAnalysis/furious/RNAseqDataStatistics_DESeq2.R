#Author : Keime, Velt
#Date : 12/2013

#########################################################################################################################
# This script performes all the comparisons described in design file, with DESeq2 library
# It uses "comp2cond_DESeq2" script if there are replicates in the experiment
# and it uses "comp2cond_DESeq2_withoutRep" script if there are no replicates
#
# comparisons are done only on samples compared, and a new normalization is done and does not match 
# that normalization values of alldata.
#
# Input: 
# - design_file: a file describing the data and desired contrastes
# - rawdata: a tab delimited file containing Ensembl gene ID and raw reads counts for each gene and each sample
# - genes_diff : path where the output tab delimited file genes_diff should be stored.
#   and log2 fold-change & p-value for each comparison describe in design file
#
# Output:
# - genes_diff.txt : A table containing annotation of all genes and raw and normalized reads counts for each gene
# - one scatterplot and one MA plot per contrast
#########################################################################################################################

RNAseqDataStatisticsDESeq2 = function(design_file, rawdata, alldata, genes_diff, allres, tlfc, tpadj, fittype){

	suppressMessages(detach(package:DESeq))
	suppressMessages(detach(package:Biobase))
	suppressMessages(library(Biobase))
	suppressMessages(library(DESeq2))

	# A table is created from the design file
	designFile=read.table(design_file, sep="\t", header=T,check.names = FALSE)	
	samplename=designFile$Sample_name
	colnames(rawdata)=samplename
	# We determined the number of column of the design file
	lengthFile=dim(designFile)[2]
	# Recovery of column corresponding to contrasts
	colContrast=designFile[,4:lengthFile]

	#Recovery of all conditions (third column of design file)
	conds = designFile[,3]
	#Recovery of number of contrasts (the first three columns correspond respectively to sampleID, samplenames and conds)
	nbrContraste=(lengthFile-3)
				
	#Design of experiment
	ExpDesign = data.frame(Condition=factor(conds,levels=levels(conds))) 
	rownames(ExpDesign)=colnames(rawdata)
			
	dds = DESeqDataSetFromMatrix(countData=rawdata,colData=ExpDesign,design=~Condition)
	dds=DESeq(dds, fitType=fittype)

	#PCA plot
	rld = rlogTransformation(dds, blind=TRUE)
	png("figures/PCA_DESeq2.png")
	suppressWarnings(plotPCA(rld, intgroup="Condition"))
	dev.off()
	
	#heatmap
	suppressMessages(library("RColorBrewer"))
	suppressMessages(library("gplots"))
	distsRL = dist(t(assay(rld)))
	mat = as.matrix(distsRL)
	rownames(mat) = colnames(mat) = with(colData(dds), Condition)
	hc=hclust(distsRL)
	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
	png("figures/heatMap_DESeq2.png")
	suppressWarnings(heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13)))
	dev.off()

	# heatmap on the 30 moste expressed genes
	library("ggplot2")
	select=order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:30]
	png("figures/heatMap_normalized_counts_DESeq2.png")
	suppressWarnings(heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6)))
	dev.off()

	png("figures/heatMap_rlogTransformation_counts_DESeq2.png")
	suppressWarnings(heatmap.2(assay(rld)[select,], col = hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6)))
	dev.off()

	tlfc = 1
	tpadj = 0.05
	
	# A loop on all contrasts column, to determine if there are replicates in the experiment
	# If there are more than one "1", there are replicates and Rep take "YES" value	
	for (i in (4:lengthFile)){

		counter = 0
		
		for (b in designFile[,i]){
			if (!is.na(b)){
				if (b == 1){
						
					counter = counter + 1 
				}
			}
		}

		# recovering of the control name (when contrast=1 in design file)
		# is used then for the "relevel" function, so that the comparison is done in the "right direction"
		# It will be convenient to make sure that Control is the first level in the treatment factor, so that the 
		# log2 fold changes are calculated as treatment over control. The function relevel achieves this:
		# > dds$Condition <- relevel( dds$Condition , "Control" )
		
		Cond1 = designFile$Condition[(designFile[,i]==1)]
		Cond1_unique = unique(Cond1)
		if (any(is.na(Cond1_unique))){
					
			Cond1_unique = Cond1_unique[-which(is.na(Cond1_unique))] #Here, it remains one condition (cond1)
					
		}
		Cond1_unique=as.character(Cond1_unique)

		Cond2 = designFile$Condition[(designFile[,i]==-1)]
		Cond2_unique = unique(Cond2)
		if (any(is.na(Cond2_unique))){
					
			Cond2_unique = Cond2_unique[-which(is.na(Cond2_unique))] #Here, it remains one condition (cond1)
					
		}

		Cond2_unique=as.character(Cond2_unique)

		compName=paste(Cond2_unique,paste("_vs_",Cond1_unique,sep=""),sep="")		

		designFileok=designFile[which(!is.na(designFile[,i])),c(0:lengthFile)]
		condsComp = designFileok[,3]
		# supress the levels which not occur in this comparison
		condsComp=factor(condsComp)
		# creation of design specific to this comparison
		ExpDesignOK = data.frame(Condition=factor(condsComp,levels=levels(condsComp)))
		rownames(ExpDesignOK)=designFileok$Sample_name
		# selection of columns to compare in rawdata
		# selection of columns to compare in rawdata
		vec = vector()
		for (v in ( 1:length(designFileok$Sample_name))) {
			int=which(colnames(rawdata)==as.character(designFileok$Sample_name[v]))
			vec[v]=int
		}
		
		# index to use to plot the scatterplot
		index_cond1_for_scatterplot = vector()
		for (index in ( 1:length(designFileok$Sample_name))) {
		if(any(grepl(as.character(Cond1_unique),designFileok[index,2]))) {
		index_cond1_for_scatterplot[index]=index
		}
		}
		index_cond1_for_scatterplot=which(!is.na(index_cond1_for_scatterplot))

		index_cond2_for_scatterplot = vector()
		for (index in ( 1:length(designFileok$Sample_name))) {
		if(any(grepl(as.character(Cond2_unique),designFileok[index,2]))) {
		index_cond2_for_scatterplot[index]=index
		}
		}
		index_cond2_for_scatterplot=which(!is.na(index_cond2_for_scatterplot))

		rawdataOK=rawdata[vec]
			
		ddsComp = DESeqDataSetFromMatrix(countData=rawdataOK,colData=ExpDesignOK,design=~Condition)
		ddsComp = DESeq(ddsComp, fitType=fittype)

		if (counter > 1){

		pdf(paste("figures/dispersion_estimation_",compName,"_DESeq2.pdf",sep=""))
		plotDispEsts(ddsComp) 
		dev.off()

		} else {

		pdf(paste("figures/dispersion_estimation_",compName,"_DESeq2_no_replicates.pdf",sep=""))
		plotDispEsts(ddsComp) 
		dev.off()

		}
		
		# if the contrast is a replicate comparison, the comparison is performed, else not
		# the comparison will performe with DESeq uniquely
		if ( i == 4){
				
			genes.res=alldata
			liste=comp2cond_DESeq2(rawdataOK, ddsComp, Cond1_unique, Cond2_unique, tlfc, tpadj, genes.res, counter, index_cond1_for_scatterplot, index_cond2_for_scatterplot)
			allres=liste$"allres"
						
		} else {
					
			genes.res=liste$"genes.res"
			liste=comp2cond_DESeq2(rawdataOK, ddsComp, Cond1_unique, Cond2_unique, tlfc, tpadj, genes.res, counter, index_cond1_for_scatterplot, index_cond2_for_scatterplot)
			allres=merge(allres, liste$"allres", by="row.names")
			rownames(allres)=allres[,1]
			allres=allres[,-1]
						
		}
				
	}
	
	genes.res=liste$"genes.res"
	write.table(genes.res, genes_diff, sep="\t", quote=F, row.names=F, dec=".")
	write.table(allres, genes_diff, sep="\t", quote=F, row.names=F, dec=".")

}





