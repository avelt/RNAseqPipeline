#
# SERE [Single-parameter quality control and sample comparison for RNA-Seq]
# 

# SERE score 
# obs.count : table of gene and sample level read counts
# TH : minimum number of reads for a gene for exclusison
sere.score <- function(obs.count, TH = 1) {
  # number of samples
	num.samp <- ncol(obs.count);
	
	# total reads across all the genes and samples
	tot.count <- sum(obs.count);
	samp.count <- colSums(obs.count);
	
	# total reads across a gene for all the samples
	row.count <- rowSums(obs.count);
	
	# filter those genes for which total reads for 
	# all the samples > TH
	idx <- which(row.count > TH);
	obs.count <- obs.count[idx, ];
	row.count <- row.count[idx];
	num.genes <- nrow(obs.count);
	
	# expected counts
	expt.count <- matrix(NA, nrow = num.genes, ncol = num.samp);
	for(gene.idx in seq(num.genes)) {
		expt.count[gene.idx, ] <- samp.count * row.count[gene.idx] / tot.count;
	}
	
	# sere score
	disp.sum <- sum((obs.count - expt.count) ^ 2 / expt.count);
	sere <- sqrt(disp.sum / (num.genes * (num.samp - 1))); 
}


# SERE dendrogram
# obs.count : table of gene and sample level read counts
# TH : minimum number of reads for a gene for exclusison
sere.dendro <- function(obs.count, TH = 1) {
	# number of samples
	num.samp <- ncol(obs.count);
	col.names <- names(obs.count);
	
	# distance matrix
	dist.mat <- matrix(NA, nrow = num.samp, ncol = num.samp);
	rownames(dist.mat) <- col.names;
	colnames(dist.mat) <- col.names;	
	
	# fill the distance matrix
	for(i in seq(num.samp)) {
		for(j in i : num.samp) {
			#print(paste('Computing SERE score for : ', col.names[i], col.names[j]));
			dist.mat[i, j] <- sere.score(obs.count[, c(i, j)], TH);
			dist.mat[j, i] <- sere.score(obs.count[, c(i, j)], TH);
		}
	}
	
	# make the dedrogram
	clust <- hclust(as.dist(dist.mat), method = "complete")
	#plot(clust)
	# done !
	return(clust);
}
