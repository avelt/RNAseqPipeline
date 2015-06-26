#Author : keime
#Date : 2011/12

###############################################################################################
#input : a data.frame with rows=genes and columns=samples, rownames=gene names, colnames=sample names
#output : a series of values and graph to explore these RNAseq data

#require : a directory called figures
#################################################################################################

RNAseqDataExploration = function(data, sizeFactor, design_file){
  #Number of libraries
  nblib = dim(data)[2]

  #Library names
  libnames = colnames(data)
 
  #Number of genes
  ################
  print("Number of genes")
  print(dim(data)[1])
  
  #% of genes expressed in each library
  print("% of genes expressed in each library")
  for (i in 1:nblib){
    #% of genes expressed in the library i
    p = round(sum(data[,i]!=0)*100/dim(data)[1], digits=2)
    print(paste(libnames[i],":",p))
  }
  
  #Density plots
  ##############
  pdf("figures/rawdata_densities.pdf")
  suppressWarnings(plot(density(log2(data[,1])), main="Rawdata densities", xlab="log2(raw read count)"))
  for (i in 2:nblib){
    	lines(density(log2(data[,i])), col=i)
  }
  legend("topright", col=1:nblib, lty=1, legend=libnames)
  invisible(dev.off())
  
  png("figures/rawdata_densities.png")
  suppressWarnings(plot(density(log2(data[,1])), main="Rawdata densities", xlab="log2(raw read count)"))
  for (i in 2:nblib){
    	lines(density(log2(data[,i])), col=i)
  }
  legend("topright", col=1:nblib, lty=1, legend=libnames)
  invisible(dev.off())
  	
  
  
  #Correlation coefficients
  #########################
  print("Pearson correlation coefficients")
  print(round(cor(data), digits=2))
  
  print("Spearman correlation coefficients")
  print(round(cor(data, method="spearman"), digits=2))
  
  # for sere.score function, it is possible to put an extra argument, TH
  # which allows to exclude genes with low reads
  print("SERE correlation coefficients")
  num.samp = dim(data)[2]
  sere.mat <- matrix(NA, nrow = num.samp, ncol = num.samp)
  for(i in seq(num.samp)) {
   for(j in i : num.samp) {
    sere.mat[i, j] <- sere.score(data[, c(i, j)]);
    sere.mat[j, i] <- sere.score(data[, c(i, j)]);
   }
  }
  colnames(sere.mat) = colnames(data)
  rownames(sere.mat) = colnames(data)
  print(round(sere.mat, digits=2))
  
  #Clustering
  ############
  
  #on raw data
  pdf("figures/rawdata_clustering_spearman.pdf")
  dist = as.dist(1-cor(data, method="spearman"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Rawdata clustering Spearman",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  pdf("figures/rawdata_clustering_pearson.pdf")
  dist = as.dist(1-cor(data, method="pearson"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings( plot(hc,main="Rawdata clustering Pearson",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  pdf("figures/rawdata_clustering_SERE.pdf")
  hc = sere.dendro(data)
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc,main="Rawdata clustering SERE"))
  invisible(dev.off())
  
  png("figures/rawdata_clustering_spearman.png")
  dist = as.dist(1-cor(data, method="spearman"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Rawdata clustering Spearman",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  png("figures/rawdata_clustering_pearson.png")
  dist = as.dist(1-cor(data, method="pearson"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc,main="Rawdata clustering Pearson",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  png("figures/rawdata_clustering_SERE.png")
  hc = sere.dendro(data)
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc,main="Rawdata clustering SERE"))
  invisible(dev.off())
 
  ###!!!
  #caution, the following clusterings are calculated with the variance estimated with the blind method
  #if there are replicates one can have a better estimation of the variance
  ###!!!
  
  #on variance stabilized data
  pdf("figures/vsddata_clustering_pearson.pdf")
  suppressMessages(library(DESeq))
  conds = factor(1:nblib)
  cds = newCountDataSet(data, conds)
  cds = estimateSizeFactors(cds) 
  cds = estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local")
  vsd = getVarianceStabilizedData(cds)
  dist = as.dist(1-cor(vsd, method="pearson"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Vsddata clustering Pearson",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  	
  pdf("figures/vsddata_clustering_euclidean.pdf")
  dists = dist(t(vsd))
  heatmap(as.matrix(dists), symm=T, main="Vsddata clustering Euclidean")
  invisible(dev.off())
  
  png("figures/vsddata_clustering_pearson.png")
  suppressMessages(library(DESeq))
  conds = factor(1:nblib)
  cds = newCountDataSet(data, conds)
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local")
  vsd = getVarianceStabilizedData(cds)
  dist = as.dist(1-cor(vsd, method="pearson"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Vsddata clustering Pearson",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  png("figures/vsddata_clustering_euclidean.png")
  dists = dist(t(vsd))
  heatmap(as.matrix(dists), symm=T, main="Vsddata clustering Euclidean")
  invisible(dev.off())
  	
  #on normalized data
  suppressMessages(library(DESeq))
  conds = factor(1:nblib)
  cds = newCountDataSet(data, conds)
  cds = estimateSizeFactors(cds)
  datanorm = t(t(data)/sizeFactors(cds))
  	
  #Density plots
  pdf("figures/normdata_densities.pdf")
  suppressWarnings(plot(density(log2(datanorm[,1])), main="Normdata densities", xlab="log2(normalized read count)"))
  for (i in 2:nblib){
  	lines(density(log2(datanorm[,i])), col=i)
  }
  legend("topright", col=1:nblib, lty=1, legend=colnames(datanorm))
  invisible(dev.off())
  
  png("figures/normdata_densities.png")
  suppressWarnings(plot(density(log2(datanorm[,1])), main="Normdata densities", xlab="log2(normalized read count)"))
  for (i in 2:nblib){
  	lines(density(log2(datanorm[,i])), col=i)
  }
  legend("topright", col=1:nblib, lty=1, legend=colnames(datanorm))
  invisible(dev.off())
  
  #Clustering
  pdf("figures/normdata_clustering_spearman.pdf")
  dist = as.dist(1-cor(datanorm, method="spearman"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Normdata clustering Spearman",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  	
  pdf("figures/normdata_clustering_pearson.pdf")
  dist = as.dist(1-cor(datanorm, method="pearson"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Normdata clustering Pearson",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  png("figures/normdata_clustering_spearman.png")
  dist = as.dist(1-cor(datanorm, method="spearman"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Normdata clustering Spearman",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  	
  png("figures/normdata_clustering_pearson.png")
  dist = as.dist(1-cor(datanorm, method="pearson"))
  hc = hclust(dist, method="average")
  if (nblib==2){
	hc=as.dendrogram(hc)
  }
  suppressWarnings(plot(hc, main="Normdata clustering Pearson",sub="",xlab="", ylab="Distance"))
  invisible(dev.off())
  
  
  #Correlation coefficients
  print("Pearson correlation coefficients on normalized data")
  print(round(cor(datanorm), digits=2))
  
  print("Spearman correlation coefficients on normalized data")
  print(round(cor(datanorm, method="spearman"), digits=2))

  
  #MDS
  #####
  
  # with replicates uniquely
   if (nblib>=4){
	  suppressMessages(library(edgeR))
	  d = DGEList(counts=data, group=1:nblib, remove.zeros=T)
	  d = calcNormFactors(d)
	  
	  #MDS plot
	  pdf("figures/mdsplot.pdf")
	  plotMDS.DGEList(d)
	  invisible(dev.off())
	  
	  png("figures/mdsplot.png")
	  plotMDS.DGEList(d)
	  invisible(dev.off())
  }
  
  #AFC on variance stabilized data
  ################################
  suppressMessages(library(ade4))
  #only if all values in vsd are > 0
  nbval = dim(vsd)[1]*dim(vsd)[2]
  if (sum(vsd>0)==nbval){ #if all values in vsd are > 0
  	pdf("figures/vsddata_afc.pdf")
  	coa = dudi.coa(vsd, scannf = F, nf = 2)
  	suppressWarnings(plot(coa$co$Comp1, coa$co$Comp2, type="n", main="Vsddata afc", xlab="Axis 1", ylab="Axis 2"))
  	text(coa$co$Comp1, coa$co$Comp2, rownames(coa$co), cex=0.8)
  	invisible(dev.off())
  	
  	png("figures/vsddata_afc.png")
  	coa = dudi.coa(vsd, scannf = F, nf = 2)
  	suppressWarnings(plot(coa$co$Comp1, coa$co$Comp2, type="n", main="Vsddata afc", xlab="Axis 1", ylab="Axis 2"))
  	text(coa$co$Comp1, coa$co$Comp2, rownames(coa$co), cex=0.8)
  	invisible(dev.off())
  	
 
	#% of explained variance
	print("% of variance explained :")
	print(round(coa$eig*100/sum(coa$eig), digits=2))
	
   }
    
  #All pairwise comparisons
  #########################

  pdf(file="figures/allxycomp.pdf")
  suppressWarnings(pairs(data, log="xy", main="allxycomp",panel=function(x,y) {abline(0,1,col="red")
  points(x,y,pch=20, cex=0.6)}))
  invisible(dev.off())

  png(file="figures/allxycomp.png")
  suppressWarnings(pairs(data, log="xy", main="allxycomp",panel=function(x,y) {abline(0,1,col="red")
  points(x,y,pch=20, cex=0.6)}))
  invisible(dev.off())

}
