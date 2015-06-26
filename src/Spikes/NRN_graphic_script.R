# Autor: Amandine Velt
# Date: 19/02/2013
##########################################################################
# Creating graphs Concentration = f(Number of reads) for all samples.
# File : path to the file htseq_combined.txt
# Dir: path to $outdir/quality_analysis/Spikes
# Length: path to the file Spikes_Length_sort.txt
# Concentration: path to the file concentration of spikes
# A graph per sample is generated. It represents the correlation between the number of reads of a spike and its concentration in the mix used.
##########################################################################

NRN = function(file, dir, length, concentration){
	
##### Data normalization of file HTSeq_combined #####
datafile = read.table(file, header = T, check.names=F, sep="\t") 
data = datafile[apply(datafile[,-1],1,sum)!=0,] # We do not take the first column (which corresponds to the id and is not numeric)
suppressMessages(library(DESeq))
nombrecol = ncol(datafile) -1 # counts the number of column of the file HTSeq (number of condition)
conds = factor(c(1:nombrecol)) # factor comprising all the conditions
cds = newCountDataSet(data[-1], conds) 
cds = estimateSizeFactors(cds)
sizefactor=sizeFactors(cds)
# if size factor = NA, replace it by 1

for (i in 1:nombrecol){
	if(is.na(sizefactor[i])) {
		sizefactor[i]=1
	}
}

datanorm = t(t(data[,-1])/sizefactor)

### Creating variables needed to create graphics Concentration = f (NRN) ###
IDgene=data[,1] # Spikes ID
tablecond = data.frame(IDgene, datanorm)
txtfile=paste(outdir,"Datanorm.txt",sep="/")
write.table(tablecond, txtfile, sep="\t", row.names = FALSE, col.names = TRUE)
LengthFile = read.table(length,header = F, check.names=F, sep=" ")
colnames(LengthFile) <- c("IDgene","Length")
TableLengthNRN = merge(LengthFile, tablecond, by = "IDgene", all = TRUE) # table containing NA for spikes which are not present
nombrecol2 = ncol(TableLengthNRN)
LengthSpikes = TableLengthNRN[,2] # length of spikes
datanormok=TableLengthNRN[ , c(3:nombrecol2)]  # data.frame containing datanorm
NRN = datanormok / LengthSpikes # create a table containing the NRN for each spike in each condition
Concentration = read.table(concentration, header=T)
nombrecol3 = ncol(datanormok) # total number of samples
ncolS1=grep("S1",names(datanormok))
ncolS2=grep("S2",names(datanormok))
NRNmix1=NRN[ , c(ncolS1)] # Select the columns of samples corresponding to Mix 1 
NRNmix2=NRN[ , c(ncolS2)] # Select the columns of samples corresponding to Mix 2 

ConcMIX1 = Concentration[,2] # Concentration of Mix1
ConcMIX2 = Concentration[,3] # Concentration of Mix2
NRNmix1ok=data.frame(NRNmix1, ConcMIX1)
NRNmix2ok=data.frame(NRNmix2, ConcMIX2)

ncol1=grep("S1",names(datafile))
ncol2=grep("S2",names(datafile))
IDGene=datafile[,1]
#Creation of a table containing sample corresponding to mix 1 only (number of reads per spike)
ReadsNumber1=datafile[ , c(ncol1)]
#Creation of a table containing sample corresponding to mix 2 only (number of reads per spike)
ReadsNumber2=datafile[ , c(ncol2)]
#Creation of a table containing Spike ID, sample corresponding to mix 2 only (number of reads per spike)
# and concentration of each spike
TableConcReadMix2=data.frame(IDGene, ReadsNumber2, ConcMIX2)
#idem
TableConcReadMix1=data.frame(IDGene, ReadsNumber1, ConcMIX1)
#Same tables but with the concentration column in ascending row
TableConcReadsMix2=TableConcReadMix2[with(TableConcReadMix2, order(ConcMIX2)), ]
TableConcReadsMix1=TableConcReadMix1[with(TableConcReadMix1, order(ConcMIX1)), ]
write.table(TableConcReadsMix1, paste(dir, "Number_Reads_Concentration_Table_Mix1.txt", sep="/"),col.names=TRUE,row.names=FALSE, quote= FALSE, sep="\t")
write.table(TableConcReadsMix2, paste(dir, "Number_Reads_Concentration_Table_Mix2.txt", sep="/"),col.names=TRUE,row.names=FALSE, quote= FALSE, sep="\t")

NbrERCC=length(tablecond$IDgene)
#Calculating the percentage of spikes present in the experiment
ERCC=(NbrERCC*100)/92
#Writing the percentage of spikes and the correlation coefficient (between number of reads and spikes concentration) 
#for each sample present within the experiment in a txt file
txt="% of spikes represented by reads (mean between all samples):"
cat(txt, file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
cat(ERCC, file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
txt2="Correlation coefficient between Concentration and NRN for each sample:"
cat(txt2, file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)

#####################GRAPHICS############################
#Creation of correlation coefficient graphe for each sample

# if the dimension of the vector is NULL, that is there is only one sample and put in one dimension
# else we take the dimension of the vector as the loop end

if (is.null(dim(NRNmix1)[2])){
	dim1=1 
} else {
	dim1=dim(NRNmix1)[2]
}

if (is.null(dim(NRNmix2)[2])){
	dim2=1
} else {
	dim2=dim(NRNmix1)[2]
}

if (dim1==1){

		pngfile3=paste(paste(dir,colnames(NRN[ncolS1]),sep="/"),"_NRN.png",sep="")
		png(pngfile3)
		dfr<-data.frame(x = NRNmix1, y = ConcMIX1)
		plot(y ~ x, xlab = "NRN", ylab = "Concentration", log = "xy", data = dfr)
		fit <- lm(formula = ConcMIX1 ~ NRNmix1)
		y1<-fit$coefficients[[2]]*0.0001
		y2<-fit$coefficients[[2]]*max(NRNmix1, na.rm=T)
		cx<-c(0.0001, max(NRNmix1, na.rm=T))
		cy<-c(y1, y2)
		lines(cx, cy)
		r=cor(NRNmix1, ConcMIX1, use = "complete.obs")
		rok=round(r,2)
		cat(colnames(NRN[ncolS1]) , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
		cat(rok , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"),sep="\n", append = TRUE)
		if(sizefactor[i]==1) {
			title(paste("Warning : Data not normalized ! / r =", rok, sep=" "))
		} else {
			title(paste("r =", rok, sep=" "))
		}
		dev.off()

} else {
	
	for (c in 1:dim(NRNmix1)[2]){ # for each column of table NRNmix1
		pngfile3=paste(paste(dir,colnames(NRNmix1)[c],sep="/"),"_NRN.png",sep="")
		png(pngfile3)
		dfr<-data.frame(x = NRNmix1[,c], y = ConcMIX1)
		plot(y ~ x, xlab = "NRN", ylab = "Concentration", log = "xy", data = dfr)
		fit <- lm(formula = ConcMIX1 ~ NRNmix1[,c])
		y1<-fit$coefficients[[2]]*0.0001
		y2<-fit$coefficients[[2]]*max(NRNmix1[,c], na.rm=T)
		cx<-c(0.0001, max(NRNmix1[,c], na.rm=T))
		cy<-c(y1, y2)
		lines(cx, cy)
		r=cor(NRNmix1[,c], ConcMIX1, use = "complete.obs")
		rok=round(r,2)
		cat(colnames(NRNmix1)[c] , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
		cat(rok , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"),sep="\n", append = TRUE)
		if(sizefactor[i]==1) {
			title(paste("Warning : Data not normalized ! / r =", rok, sep=" "))
		} else {
			title(paste("r =", rok, sep=" "))
		}
		dev.off()
	}
	
}

if (dim2==1){

	pngfile4=paste(paste(dir,colnames(NRN[ncolS2]),sep="/"),"_NRN.png",sep="")
	png(pngfile4)
	dfr<-data.frame(x = NRNmix2, y = ConcMIX2)
	plot(y ~ x, xlab = "NRN", ylab = "Concentration", log = "xy", data = dfr)
	fit <- lm(formula = ConcMIX2 ~ NRNmix2)
	y1<-fit$coefficients[[2]]*0.0001
	y2<-fit$coefficients[[2]]*max(NRNmix2, na.rm=T)
	cx<-c(0.0001, max(NRNmix2, na.rm=T))
	cy<-c(y1, y2)
	lines(cx, cy)
	r=cor(NRNmix2, ConcMIX2, use = "complete.obs")
	rok=round(r,2)
	cat(colnames(NRN[ncolS2]) , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
	cat(rok , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
	if(sizefactor[i]==1) {
		title(paste("Warning : Data not normalized ! / r =", rok, sep=" "))
	} else {
		title(paste("r =", rok, sep=" "))
	}
	dev.off()

} else {
	
	for (c in 1:dim(NRNmix2)[2]){
		pngfile4=paste(paste(dir,colnames(NRNmix2)[c],sep="/"),"_NRN.png",sep="")
		png(pngfile4)
		dfr<-data.frame(x = NRNmix2[,c], y = ConcMIX2)
		plot(y ~ x, xlab = "NRN", ylab = "Concentration", log = "xy", data = dfr)
		fit <- lm(formula = ConcMIX2 ~ NRNmix2[,c])
		y1<-fit$coefficients[[2]]*0.0001
		y2<-fit$coefficients[[2]]*max(NRNmix2[,c], na.rm=T)
		cx<-c(0.0001, max(NRNmix2[,c], na.rm=T))
		cy<-c(y1, y2)
		lines(cx, cy)
		r=cor(NRNmix2[,c], ConcMIX2, use = "complete.obs")
		rok=round(r,2)
		cat(colnames(NRNmix2)[c] , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
		cat(rok , file =paste(dir,"coefficient_correlation_percent_spikes.txt",sep="/"), sep="\n", append = TRUE)
		if(sizefactor[i]==1) {
			title(paste("Warning : Data not normalized ! / r =", rok, sep=" "))
		} else {
			title(paste("r =", rok, sep=" "))
		}
		dev.off()
	}
	
}


}




