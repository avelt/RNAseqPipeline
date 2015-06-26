#Author : Keime
#Date : 02/2013

########################################################################################################
#input : - a data.frame with rows=genes and columns=samples, rownames=gene names, colnames=sample names
# - dist_method : The distance measure to be used. This must be one of "euclidean", "maximum", 
# "pearson", "correlation", "spearman" or "kendall". 
# - hclust_method : There are several alternative clustering methods.
# We use average linkage method single (others possible : linkage, complete linkage...)
#
#output : a graphic of clustering of all samples
#
########################################################################################################

RNAseqClustering = function(rawdata, dist_method, hclust_method){

nblib = dim(rawdata)[2]
libnames = colnames(rawdata)
suppressMessages(library(DESeq))
conds = factor(1:nblib)
cds = newCountDataSet(rawdata, conds)
cds = estimateSizeFactors(cds)
datanorm = t(t(rawdata)/sizeFactors(cds))
dist = as.dist(1-cor(datanorm, method=dist_method))
hc = hclust(dist, method=hclust_method)

if (length(conds)==2){
	hc=as.dendrogram(hc)
}

plot(hc, main="",sub="",xlab="", ylab="Distance")

}
