# Autor: Amandine Velt
# Date: 19/02/2013
#########################################################################
# Creating graph Expected Fold Change=f(Observed Fold Change) for all samples.
# File : path to the file htseq_combined.txt
# Outdir: path to $outdir/quality_analysis/Spikes
# foldexpected: path to the file Expected_FoldChange.txt
# A graph representing the expected fold change between the mix of spikes as a function of observed fold change in the RNA-seq experiment is plotted. One graph represents all the samples. There are 92 points, corresponding to 92 spikes.
##########################################################################


FoldChange = function(file, outdir, foldexpected){
 
NumberReadsMix=read.table(file, header=T, check.names=F, sep="\t") 

ncolS1=grep("S1",names(NumberReadsMix)) #selects the columns containing S1 in their ID
ncolS2=grep("S2",names(NumberReadsMix)) #selects the columns containing S2 in their ID
ReadsNumberS1=NumberReadsMix[ , c(ncolS1)] #Select the columns of samples corresponding to Mix 1.
ReadsNumberS2=NumberReadsMix[ , c(ncolS2)] # Select the columns of samples corresponding to Mix 2. 

if ((length(ncolS1)==1) & (length(ncolS2)==1)){

	ObservedFoldChange = ReadsNumberS1/ReadsNumberS2

} else {

	ReadsSumMix1=rowSums(ReadsNumberS1) # Sum of lines (sum of the number of reads of all samples corresponding to Mix 1)
	ReadsMeanMix1=ReadsSumMix1 /2 # and calculating the mean.
	ReadsSumMix2=rowSums(ReadsNumberS2)
	ReadsMeanMix2=ReadsSumMix2 /2
	ObservedFoldChange = ReadsMeanMix1 / ReadsMeanMix2

}

IDgene=NumberReadsMix[,1]
ObservedFoldchange = data.frame(IDgene, ObservedFoldChange)
ExpectedFoldChange=read.table(foldexpected, header=T)

TableFoldChange = merge(ObservedFoldchange, ExpectedFoldChange, by = "IDgene", all = TRUE)

FoldChangeO=TableFoldChange[,2] # Recovery of observed fold-change
FoldChangeE=TableFoldChange[,3] # Recovery of expected fold-change

txtfile=paste(outdir,"Table_Fold_Change_Spikes.txt",sep="/")

write.table(TableFoldChange, txtfile, sep="\t", row.names = FALSE, col.names = TRUE) 

# Creation of graphic of Fold Change
pngfile4=paste(outdir,"foldChange_spikes.png",sep="/")
png(pngfile4)
colors = rep(1,92)
if ((length(ncolS1)==1) & (length(ncolS2)==1)){
colors[ReadsNumberS1<50] = 10
} else {
colors[ReadsMeanMix1<50] = 10 # red dots representing the ERCC with less than 50 reads
}
plot( x = FoldChangeO,
y = FoldChangeE,
xlab = "Observed Fold Change",
ylab = "Expected Fold Change",
xlim = (c(0,5)),
ylim = (c(0,5)),
pch = 16,
col = colors)
legend("topleft", legend = c("Reads Number < 50", "Reads Number > 50"), col = c("red", "black"), pch = 16)
abline(0,1)
dev.off()
}


