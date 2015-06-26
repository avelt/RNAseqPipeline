#Author : Velt
#Date : 02/2013

#########################################################################################################################
# This script performes all the comparisons described in design file, with DESeq library
# It uses "comp2cond.R" script if there are replicates in the experiment
# and it uses "comp2cond_var2samples.R" script if there are no replicates
#
# Input: 
# - design_file: a file describing the data and desired contrastes
# - rawdata: a tab delimited file containing Ensembl gene ID and raw reads counts for each gene and each sample
# - genes_diff : path where the output tab delimited file genes_diff should be stored.
#   and log2 fold-change & p-value for each comparison describe in design file
# - allres_file : path where the output tab delimited file allres should be stored.
#
# Output:
# - genes_diff.txt : A table containing annotation of all genes and raw and normalized reads counts for each gene
# - allres.txt : A table containing Ensembl gene ID and one column per sample
#   The boxes values can be 0, 1 or -1
#   0 means that gene is not differentially expressed
#	1 means that gene is not overexpressed in condition B regarding to condition A
#	-1 means that gene is not underexpressed in condition B regarding to condition A
# - one scatterplot per contrast
#########################################################################################################################

RNAseqDataStatistics = function(design_file, rawdata, genes_diff, allres_file, tlfc, tpadj){
	
	# A table is created from the design file
	designFile=read.table(design_file, sep="\t", header=T,check.names = FALSE)							
	# We determined the number of column of the design file
	lengthFile=dim(designFile)[2]
	# Recovery of column corresponding to contrasts
	colContrast=designFile[,4:lengthFile]
	# Default value of replicate
	Rep="NO"
	
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

		if (counter > 1){

		Rep="YES"
				
		}

	}

		
	counter=0
	#Recovery of all conditions (third column of design file)
	conds = designFile[,3]
	conds=as.factor(conds)
	#Recovery of number of genes in rawdata to create allres table
	nbgenes = dim(rawdata)[1]
	#Recovery of number of contrasts (the first three columns correspond respectively to sampleID, samplenames and conds)
	nbrContraste=(lengthFile-3)
	#Creation of allres table and attribution of rownames
	allres = matrix(data = 0, nrow = nbgenes, ncol=nbrContraste)
	allresOk = as.data.frame(allres)
	rownames(allresOk) = rownames(rawdata)
	
	#If there are replicates in the experiment
	if (Rep=="YES"){
		
		#Create a table with rownames=ensembl gene ID, colnames=conds
		#and the number of reads in each gene and each condition
		cds = newCountDataSet(rawdata, conds)
		#Estimate the size factors from the count data
		cds = estimateSizeFactors( cds )
		#Estimate the dispersion on all samples
		cds = estimateDispersions( cds )
		
		#Plot the variance on all samples
		plotDispEsts = function( cds ){
			suppressWarnings(plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)$perGeneDispEsts, pch = '.', log="xy"))
			xg = 10^seq( -.5, 5, length.out=300 )
			lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
		}
			
		pdf("figures/variance.pdf")
		plotDispEsts( cds )
		invisible(dev.off()) 
		
		#Here, pairwise comparisons begin.
		#Contrasts are treated one to one
		for (i in (4:lengthFile)){
			
			#Counter allow to determine in each column of allres, the results should be stored.
			#Column 1 for the first contrast, column 2 for the second ...
			counter=counter+1
			#Recovery of first condition (selection of condition(s) where contrast = 1)
			Cond1 = designFile$Condition[(designFile[,i]==1)]
			#With replicates, several conditions correspond to value 1
			#These conditions have the same name (because they are replicates)
			#So, one condition is kept 
			Cond1_unique = unique(Cond1)
			#Certain conditions have value "NA" (because they are not compared)
			#So, these conditions are suppressed
			
			if (any(is.na(Cond1_unique))){
			
				Cond1_unique = Cond1_unique[-which(is.na(Cond1_unique))] #Here, it remains one condition (cond1)
			
			}
			
			
			#Recovery of second condition (selection of condition(s) where contrast = -1)
			Cond2 = designFile$Condition[(designFile[,i]==-1)]
			Cond2_unique = unique(Cond2)
			
			if (any(is.na(Cond2_unique))){
			
				Cond2_unique = Cond2_unique[-which(is.na(Cond2_unique))] #Here, it remains one condition (cond2)
			
			}
			
			#When i = 4, it means that it is the first comparison (first column of contrast)
			#So, genes.res=alldata (the table containing no comparison)
			if ( i == 4){
						
				genes.res=alldata
				cond1=Cond1_unique
				cond2=Cond2_unique
				nb=counter
				allres=allresOk
				#generate baseMean for condition A and B, log2 fold-change and p-value between these two conditions
				#this columns are added to alldata
				listres = comp2cond(cds, as.character(cond1), as.character(cond2), allres, nb, genes.res, tlfc, tpadj)
						
			} else {
					
				cond1=Cond1_unique
				cond2=Cond2_unique
				nb=counter
				allres=listres$allres
				#we use the last modified table, containing the columns of all perfomed comparisons
				#and at each comparison, news columns are added
				genes.res=listres$genes.res
				listres = comp2cond(cds, as.character(cond1), as.character(cond2), allres, nb, genes.res, tlfc, tpadj)
							
				}
						
			}
			
			#After the last comparison, allres and genes.res are update.
			#Each table contains the results of all comparison
			allres=listres$allres
			genes.res=listres$genes.res
			
			#Each table of results is write in a file
			write.table(genes.res, genes_diff, sep="\t", quote=F, row.names=F, dec=",")
			write.table(allres, allres_file, sep="\t", quote=F, row.names=F, dec=",")
	
	#Else there are no replicates in the experiment		
	} else {
		
		for (i in (4:lengthFile)){
			
			cds = newCountDataSet(rawdata, conds)
			#Estimate the size factors from the count data
			#method="blind" and sharingMode="fit-only" are used because there are no replicates
			cds = estimateSizeFactors( cds , method="blind", sharingMode="fit-only" )
			#Estimate dispersion
			cds = estimateDispersions( cds , method="blind", sharingMode="fit-only" )
	
			plotDispEsts = function( cds ){
				suppressWarnings(plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)$perGeneDispEsts, pch = '.', log="xy"))
				xg = 10^seq( -.5, 5, length.out=300 )
				lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
			}
			
			pdf("figures/variance.pdf")
			plotDispEsts( cds )
			invisible(dev.off()) 
			
			#When there are no replicates, the dispersion is estimated between the two compared conditions
			#And not between all the samples, as previously.
			#So, the script comp2cond.var2samples takes such as argument the number corresponding to the sample
			#For example, if Cond1 corresponds to the third sample in the samplename column, nb1=3 etc ...
			#And this allow to select the two good samples in cds to estimate dispersion between their
			
			cds = newCountDataSet(rawdata, conds)
			#Estimate the size factors from the count data
			#method="blind" and sharingMode="fit-only" are used because there are no replicates
			cds = estimateSizeFactors( cds , method="blind", sharingMode="fit-only" )
			
			#The rest of code performed as above, see precedent commments
			counter=counter+1
			Cond1 = designFile$Condition[(designFile[,i]==1)]
			nb1=which(designFile[,i]==1)
			Cond1_unique = unique(Cond1)
			
			if (any(is.na(Cond1_unique))){
			
				Cond1_unique = Cond1_unique[-which(is.na(Cond1_unique))]
			
			}
			
			Cond2 = designFile$Condition[(designFile[,i]==-1)]
			nb2=which(designFile[,i]==-1)
			Cond2_unique = unique(Cond2)
			
			if (any(is.na(Cond2_unique))){
			
				Cond2_unique = Cond2_unique[-which(is.na(Cond2_unique))]
			
			}
			
					
			if ( i == 4){
						
				genes.res=alldata
				cond1=Cond1_unique
				cond2=Cond2_unique
				nb=counter
				allres=allresOk
				listres = comp2cond.var2samples(cds, nb1, nb2, as.character(cond1), as.character(cond2), allres, nb, genes.res, tlfc, tpadj)
						
			} else {
					
				cond1=Cond1_unique
				cond2=Cond2_unique
				nb=counter
				genes.res=listres$genes.res
				allres=listres$allres
				listres = comp2cond.var2samples(cds, nb1, nb2, as.character(cond1), as.character(cond2), allres, nb, genes.res, tlfc, tpadj)
						
			}
						
		}
			
		allres=listres$allres
		genes.res=listres$genes.res

		write.table(genes.res, genes_diff, sep="\t", quote=F, row.names=F, dec=".")
		write.table(allres, allres_file, sep="\t", quote=F, row.names=rownames(allres), dec=".")
							
	}
	
}
