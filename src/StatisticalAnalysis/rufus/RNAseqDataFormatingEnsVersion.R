#Author : keime
#Date : 2013/02


###############################################################################################
#input : a file with rows=genes and columns=samples (except 1st column=ensembl gene id and 1st line=library names), tab delimited
#output : a data.frame with the following columns :
#ensembl gene id
#raw read counts for each library (one column per library)
#normalized data for each library (one column per library)
#normalized data divided by gene length for each library (one column per library)
#Gene name
#Description

#require : biomart library
#################################################################################################

#arguments :
#file : path/to/file with rows=genes and columns=samples (except 1st column=ensembl gene id and 1st line=library names), tab delimited
#bmdataset :  the name of the dataset to use in Ensembl to search for corresponding annotations. For a list of all available dataset use useMart('ensembl');listDatasets(mart)
#samplenames : name of the samples (in the same order as in the file file)

RNAseqDataFormatingEnsVersion = function(file, bmdataset, host, samplenames, geneBiotypeFile, gtffile){
	
	#read data file
	datafile = read.table(file, row.names=1, h=T, check.names=F, sep="\t")
	
	#sample names
	colnames(datafile) = samplenames
	
	#suppress genes not expressed in all samples
	data = datafile[apply(datafile,1,sum)!=0,]
	
	#number of libraries
	nblib= dim(data)[2]
	
	#libraries names
	libnames=colnames(data)
	
	######################################
	######################################
	# GENE ANNOTATION
	######################################
	
	#Add annotations and calculate RPKM
	suppressMessages(library(biomaRt))
	ensembl=useMart("ENSEMBL_MART_ENSEMBL", host=host, dataset=bmdataset)
	# the 'external_gene_id' attribute has been changed to 'external_gene_name' since the 76 release so after aug2014, we use the external_gene_name attribute.
	if(length(grep("2016",host))==0 && length(grep("2015",host))==0 && length(grep("aug2014",host))==0 && length(grep("oct2014",host))==0 && length(grep("dec2014",host))==0){
		annotation1 = getBM(attributes=c("ensembl_gene_id","external_gene_id","description","chromosome_name","start_position","end_position", "gene_biotype", "ensembl_transcript_id"), 
				filters="ensembl_gene_id", values=rownames(data), mart=ensembl)
		#because all the annotations are not always found in a first step 
		not = rownames(data)[!rownames(data) %in% unique(annotation1$ensembl_gene_id)]
		if (length(not) !=0){
			annotationnot = getBM(attributes=c("ensembl_gene_id","external_gene_id","description","chromosome_name","start_position","end_position", "gene_biotype", "ensembl_transcript_id"), 
					filters="ensembl_gene_id", values=not, mart=ensembl)
			
			#concatenate all results
			annotation = rbind(annotation1, annotationnot)		
		} else {
			annotation = annotation1
		}

	} else {

		annotation1 = getBM(attributes=c("ensembl_gene_id","external_gene_name","description","chromosome_name","start_position","end_position", "gene_biotype", "ensembl_transcript_id"), 
				filters="ensembl_gene_id", values=rownames(data), mart=ensembl)
		#because all the annotations are not always found in a first step 
		not = rownames(data)[!rownames(data) %in% unique(annotation1$ensembl_gene_id)]
		if (length(not) !=0){
			annotationnot = getBM(attributes=c("ensembl_gene_id","external_gene_name","description","chromosome_name","start_position","end_position", "gene_biotype", "ensembl_transcript_id"), 
					filters="ensembl_gene_id", values=not, mart=ensembl)
			
			#concatenate all results
			annotation = rbind(annotation1, annotationnot)		
		} else {
			annotation = annotation1
		}

	}

	######################################
	######################################
	# GO ANNOTATION
	######################################
	
	#Search for all GO annotations
	go = getBM(attributes=c("ensembl_gene_id","name_1006","namespace_1003"), 
				filters="ensembl_gene_id", values=rownames(data), mart=ensembl)
	
	#matrix with GO annotations of all genes (4 columns : gene id, then the 3 ontologies)
	tabannot = matrix(nrow=length(unique(go$ensembl_gene_id)), ncol=4) 
	#the 3 different ontologies
	ontologies = c("biological_process","molecular_function","cellular_component") 
	colnames(tabannot) = c("Ensembl gene id", ontologies)
	#number of the gene under consideration
	gnum=0

	for (g in unique(go$ensembl_gene_id)){ #for all genes
		gnum = gnum+1 #line in which all information about this gene will be displayed
		tabannot[gnum,1] = g #add the ID of the gene in this line, 1st column
		for (o in ontologies){ #for the 3 different ontologies
			allannot = go[go$ensembl_gene_id==g & go$namespace_1003==o, 2] #all the annotations corresponding to this gene and this ontology
			allannotsameline = "" #character that will contain all these annotations in one line
			for (a in allannot){allannotsameline=paste(allannotsameline, a, sep="; ")} #put all the annotations in the same line, separated by "; "
			allannotsameline = substr(allannotsameline, 3, nchar(allannotsameline)) #remove the first "; "
			tabannot[gnum,o] = allannotsameline #put the annotations in tabannot (in line gnum and for the ontology o)
		}
	}
	
	#Suppress all the terms "biological process", "molecular function" or "cellular component"
	tabannot = gsub("biological_process; ","", tabannot)
	tabannot = gsub("biological_process","", tabannot)
	tabannot = gsub("molecular_function; ","", tabannot)
	tabannot = gsub("molecular_function","", tabannot)
	tabannot = gsub("cellular_component; ","", tabannot)
	tabannot = gsub("cellular_component","", tabannot)
	colnames(tabannot)[1]="ensembl_gene_id"
	#Add GO annotations to the data
	annotationgene = unique(annotation[,1:7])
	# merge of genes annotations and go annotations
	tmp = merge(annotationgene, tabannot, by="ensembl_gene_id", all.x=T)
	colnb = dim(tmp)[2]
	colnames(tmp)[(colnb-6):colnb] = c("chromosome name","start gene position","end gene position","Gene biotype","GO:biological process", "GO:molecular function", "GO:cellular component")

	#Replace NA by "" in all GO annotations
	tmp[,(colnb-2):colnb][is.na(tmp[,(colnb-2):colnb])] = ""

	data.goannotated.annotationgene = tmp
	
	######################################
	######################################
	# CALCULATION OF GENES SIZE
	######################################
	# explanation :
	# we do not want to take into account the length of lncRNAs transcripts in computing the size of the genes
	# thus the calculation of size of genes takes place in two steps
	# Indeed, the problem is that in removing these lncRNA the size of some genes (~ 500) can not be calculated
	# because some genes have only lncRNA transcripts
	# FIRST STEP : creates a vector with all the genes length (with lncRNA)
	# SECOND STEP : creates a vector with the genes length without lncRNA length
	# generation of a vector containing the length calculated in step 2 (without lncRNA), 
	# and the length of genes which code only for lncRNA
	
	# calculating the length of the genes
	suppressMessages(library(GenomicRanges))
	suppressMessages(library(rtracklayer))
	GTFfile = gtffile
	# reading of gtf file, containing ALL genes and transcripts
	GTF=import.gff(GTFfile, format="gtf", asRangedData=F, feature.type="exon")
	
	# removes retained_intron and sense_intron, because they are long non-coding transcripts
	# they must not be taken into account for the calculation of the genes size
	GTFtmp=GTF[which(GTF$source!="retained_intron")]
	GTFtmp=GTFtmp[which(GTFtmp$source!="sense_intronic")]
	# FIRST STEP
	# creates a list, with one vector per gene ID, each gene ID containing all the exon intervals.
	# the exons that overlap are grouped in a set of exons (one interval)
	grl=reduce(split(GTF, elementMetadata(GTF)$gene_id))
	reducedGTF=unlist(grl, use.names=T)
	# create a vector containing all gene ID
	elementMetadata(reducedGTF)$gene_id=rep(names(grl), elementLengths(grl))
	# create a vector containing all corresponding length for each gene ID
	elementMetadata(reducedGTF)$widths=width(reducedGTF)
	# function which allows to do sum of all exons length for one gene ID.
	calc_length=function(x) {
		sum(elementMetadata(x)$widths)
	}
	# using of calc_length function to calculate the gene length for each gene ID
	
	# certain genes reflect only long non-coding transcripts, and we're not getting the size.
	# we get here their size and we concatenate the two files
	# geneLength1 contains all length, lnc + the others
	geneLength1=as.data.frame(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
	colnames(geneLength1)=c("Gene_length")
	

	# SECOND STEP
	# creates a list, with one vector per gene ID, each gene ID containing all the exon intervals.
	# the exons that overlap are grouped in a set of exons (one interval)
	grltmp=reduce(split(GTFtmp, elementMetadata(GTFtmp)$gene_id))
	reducedGTFtmp=unlist(grltmp, use.names=T)
	# create a vector containing all gene ID
	elementMetadata(reducedGTFtmp)$gene_id=rep(names(grltmp), elementLengths(grltmp))
	# create a vector containing all corresponding length for each gene ID
	elementMetadata(reducedGTFtmp)$widths=width(reducedGTFtmp)
	# geneLength2 contains all length except length of lncRNA 
	geneLength2=as.data.frame(sapply(split(reducedGTFtmp, elementMetadata(reducedGTFtmp)$gene_id), calc_length))
	colnames(geneLength2)=c("Gene_length")
	
	# recovery of geneLength1 indices whose geneID is common with geneLength2
	# with : which(rownames(geneLength1) %in% rownames(geneLength2))
	# and supression of these indices in geneLength1 to keep only gene length which code uniquely for lncRNA
	#geneLength1$geneID=rownames(geneLength1)
	#geneLength1tmp=geneLength1[-which(rownames(geneLength1) %in% rownames(geneLength2)),] 
	geneLength2$geneID=rownames(geneLength2)
	geneLength1$geneID=rownames(geneLength1)
	geneLength=rbind(geneLength2, geneLength1[-which(rownames(geneLength1) %in% rownames(geneLength2)),])
	
	geneLengthTMP=as.data.frame(geneLength[,-2])
	rownames(geneLengthTMP)=geneLength[,2]
	geneLength=geneLengthTMP

	
	#data with gene length
	data=data[order(rownames(data)),]
	datalen = merge(data, geneLength, by="row.names") 
	colnames(datalen) = c("Ensembl_gene_id",colnames(data), "Gene_length")
	
	#data with annotations and gene length
	dataannot = merge(datalen, data.goannotated.annotationgene, by.x="Ensembl_gene_id", by.y="ensembl_gene_id")
	
	#to keep only the first part of the gene description (before [)
	tmpdesc = strsplit(as.character(dataannot$description),"[", fixed=T)
	f = function(l){
		if (length(l)>=1){
			return(l[[1]])
		} else {
			return("")
		}
	}
	tmpdescok = unlist(lapply(tmpdesc, f))
	dataannot$description = tmpdescok
	
	######################################
	######################################
	# RPKN CALCULATION
	######################################
	#nbcol = dim(dataannot)[2] #nb of column in the data.frame
	#nbalign = apply(dataannot[-c(1,nbcol,nbcol-1,nbcol-2)], 2, sum) #nb of reads aligned and assigned to a gene
	#rpkm = t  (t(dataannot[,-c(1,nbcol,nbcol-1,nbcol-2)] / (dataannot[,nbcol-2]/1000) ) / (nbalign/1000000))
	
	#normalized data calculation
	nbcol = dim(dataannot)[2] #nb of column in the data.frame
	suppressMessages(library(DESeq))
	conds = factor(1:nblib)
	cds = newCountDataSet(dataannot[,-c(1,nbcol,nbcol-1,nbcol-2,nbcol-3, nbcol-4, nbcol-5, nbcol-6,nbcol-7,nbcol-8,nbcol-9)], conds)
	cds = estimateSizeFactors(cds)
	datanorm = t(t(dataannot[,-c(1,nbcol,nbcol-1,nbcol-2,nbcol-3, nbcol-4, nbcol-5, nbcol-6,nbcol-7,nbcol-8,nbcol-9)])/sizeFactors(cds))

	#normalized data adjusted for gene length (normalized data / gene length)
	rpkn = datanorm / (as.vector(dataannot[,nbcol-9]/1000 ))
	
	#data + annotations + rpkn
	dataall = data.frame(dataannot[,-c(nbcol,nbcol-1,nbcol-2,nbcol-3, nbcol-4, nbcol-5, nbcol-6,nbcol-7,nbcol-8,nbcol-9)] , datanorm, rpkn, dataannot[,c(nbcol-9,nbcol-8,nbcol-7,nbcol-6,nbcol-5, nbcol-4, nbcol-3,nbcol-2,nbcol-1,nbcol)]  )	
	#renames columns
	colnames(dataall) = c("Ensembl gene id", paste(libnames,"(raw read counts)"), paste(libnames,"(normalized)"), paste(libnames,"(normalized and divided by gene length in kb)"), "Gene length", "Gene name", "Description","Chromosome name", "Start gene position", "End gene position", "Gene biotype", "GO:biological process", "GO:molecular function", "GO:cellular component")
	
	######################################
	######################################
	# STATISTICS ON GENE'S BIOTYPES
	######################################
	
	nbcol=dim(dataall)[2]
	proteinCoding=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="protein_coding"])/dim(dataall)[1])*100, digits=3)
	lincRNA=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="lincRNA"])/dim(dataall)[1])*100, digits=3)
	snoRNA=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="snoRNA"])/dim(dataall)[1])*100, digits=3)
	antisense=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="antisense"])/dim(dataall)[1])*100, digits=3)
	miRNA=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="miRNA"])/dim(dataall)[1])*100, digits=3)
	pseudogene=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="pseudogene"])/dim(dataall)[1])*100, digits=3)
	snRNA=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="snRNA"])/dim(dataall)[1])*100, digits=3)
	rRNA=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="rRNA"])/dim(dataall)[1])*100, digits=3)
	process=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="processed_transcript"])/dim(dataall)[1])*100, digits=3)
	misc=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="misc_RNA"])/dim(dataall)[1])*100, digits=3)
	senseIntronic=round((length(dataall[,nbcol-3][dataall[,nbcol-3]=="sense_intronic"])/dim(dataall)[1])*100, digits=3)

	cat(paste(paste("protein_coding :", proteinCoding, sep=" "), "%", sep=" "),file=geneBiotypeFile,sep="\n")
	cat(paste(paste("pseudogene :", pseudogene, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("lincRNA :", lincRNA, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("antisense :", antisense, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("processed_transcript :", process, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("snoRNA :", snoRNA, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("miRNA :", miRNA, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("snRNA :", snRNA, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("rRNA :", rRNA, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("misc_RNA :", misc, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")
	cat(paste(paste("sense_intronic :", senseIntronic, sep=" "), "%", sep=" "),file=geneBiotypeFile,append=TRUE,sep="\n")

	return(dataall)
	
}
