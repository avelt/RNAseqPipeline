#! /BIN/sh
############################################################################################################
#This script generates a wig file calculated on uniquely aligned reads and normalized
#usage : sh generateWig.sh file.sam normfactor trackline > file.wig.gz
#file.sam : a SAM file containing alignement results
#normfactor : a normalization factor to normalize coverage data
#trackline : first line of WIG file that will be generated
#file.wig : a wig file containing normalized coverage data calculated on uniquely aligned reads, gzip compressed
############################################################################################################
# exemple track type :
# nameSample=$( grep ${name} ${design_file} | cut -f2 )
# 'track type=bedGraph name="${name}" description="${nameSample}" visibility=full alwaysZero=on windowingFunction=maximum'
#arguments

samfile=${1}
normfactor=${2}
trackline=${3}
paramfile=${4}
outdir=${5}
# check of variables
echo "Variables used :"
echo "Sam file : ${samfile}"
echo "Normalization factor : ${normfactor}"
echo "Track line : ${trackline}"

source ${paramfile}
cd $outdir
#name of file (last file of the path, without extension and without "alignment")
name=`echo "${samfile}" | perl -e '$in=<>;@split=split(/\//,$in);@split2=split(/\.sam/,$split[$#split]); print substr($split2[0],0,length($split2[0])-10);'`

echo ${name}


# if the norm factor is equal to 1 we specify that the wig file will not be normalized.
if [ ${normfactor} = 1 ]
then

	wigFileName="${outdir}/Alignment/${name}/${name}_uniq_coverage_nonNormalized.wig"
	tdfFileName="${outdir}/Alignment/${name}/${name}_uniq_coverage_nonNormalized.tdf"
	echo "Normalization factor equal : ${normfactor} => wig file not normalized."

else

	wigFileName="${outdir}/Alignment/${name}/${name}_uniq_normcoverage.wig"
	tdfFileName="${outdir}/Alignment/${name}/${name}_uniq_normcoverage.tdf"

fi

echo "#####################"
echo " wig file name : ${wigFileName}"

#write track line into a file
touch ${outdir}/Alignment/${name}/tmptrack_${name}.txt
echo ${trackline} > ${outdir}/Alignment/${name}/tmptrack_${name}.txt

#recover of reference fasta file to generate bam file from sam file with -T option
#select only uniquely aligned reads
echo "`date`: ${name} - select only uniquely aligned reads in sam file "
grep "NH:i:1	" ${samfile} > ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam
grep "	NH:i:1$" ${samfile} >> ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam
wait
wait

if [ -f ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam ]
then

	# recovery of the header of the raw bam file
	# order to concatenate it with the sam containing unique reads
	# necessary to convert the new sam file in new bam file
	${BIN}/${SAMTOOLS_VERSION}/samtools view -H ${outdir}/Alignment/${name}/${name}_alignment.bam > ${outdir}/Alignment/${name}/tmp_header.sam
	# recovering the list of chromosomes present in the header
	listChr=$( cut -f2 ${outdir}/Alignment/${name}/tmp_header.sam | cut -d":" -f2 | grep '\(chr\|ERCC\)' )
	# deletion of unnecessary lines
	sed -i '/ID:TopHat/d' ${outdir}/Alignment/${name}/tmp_header.sam
	# this loop removes of the header, the chromosomes not present in the new file
	# otherwise samtools returns an error when converting sam file into bam file
	for chr in in ${listChr}
	do

		if grep -q ${chr} ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam
		then

			echo "Found" > /dev/null
			
		else
		
			sed -i '/'"${chr}"'/d' ${outdir}/Alignment/${name}/tmp_header.sam
		
		fi

	done
	wait
	wait
	# concatenation of good header with the new sam file (containing unique reads)
	cat ${outdir}/Alignment/${name}/tmp_header.sam ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam > ${outdir}/Alignment/${name}/tmpsamuniq_${name}_header.sam
	mv ${outdir}/Alignment/${name}/tmpsamuniq_${name}_header.sam ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam
	rm ${outdir}/Alignment/${name}/tmp_header.sam
	wait


	# convert sam file (unique reads) to bam (unique reads) because it's necessary for igvtools
	echo "`date`: ${name} - convert sam file containing uniquely aligned reads to bam file 
	${BIN}/${SAMTOOLS_VERSION}/samtools view -bS ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam > ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam"
	${BIN}/${SAMTOOLS_VERSION}/samtools view -bS ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam > ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam
	wait

	# sort of bam file containing only uniquely aligned reads
	echo "`date`: ${name} - sort bam file
	${BIN}/${SAMTOOLS_VERSION}/samtools sort ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam ${outdir}/Alignment/${name}/tmpbamuniq_${name}_sorted"
	${BIN}/${SAMTOOLS_VERSION}/samtools sort ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam ${outdir}/Alignment/${name}/tmpbamuniq_${name}_sorted
	mv ${outdir}/Alignment/${name}/tmpbamuniq_${name}_sorted.bam ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam
	wait
	# index bam file to generate bam.bai
	echo "`date`: ${name} - index bam file "
	${BIN}/${SAMTOOLS_VERSION}/samtools index ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam
	wait
	#just for verification
	#wc -l tmpsamuniq."$name"
	#generate wig file on uniquely aligned reads

	cd ${BIN}/${IGV_VERSION}
	
	if [ ${GENOME_VERSION} == "zv9" ]
	then
	
		GENOME_VERSION="danRer7"
	fi
		
	if [ -f ${BIN}/${IGV_VERSION}/genomes/${GENOME_VERSION}.genome ] 
	then

		#generation of wig file from bam file containing only uniquely aligned reads
		./igvtools count -w ${WINDOW_SIZE} ${outdir}/Alignment/${name}/tmpbamuniq_${name}.bam ${outdir}/Alignment/${name}/tmpwiguniq_${name}.wig ${BIN}/${IGV_VERSION}/genomes/${GENOME_VERSION}.genome
		wait
		sed -i '/#Columns: Pos, ComBINed Strands/d' ${outdir}/Alignment/${name}/tmpwiguniq_${name}.wig
		#normalize wig file
		wait
			
		echo "`date`: ${name} - perl ${SCRIPTSDIR}/src/Utilities/normalizeWig.pl ${outdir}/Alignment/${name}/tmpwiguniq.${name}.wig ${normfactor} "`cat ${outdir}/Alignment/${name}/tmptrack_${name}.txt`" "

		perl ${SCRIPTSDIR}/src/Utilities/normalizeWig.pl ${outdir}/Alignment/${name}/tmpwiguniq_${name}.wig ${normfactor} "`cat ${outdir}/Alignment/${name}/tmptrack_${name}.txt`" > ${wigFileName}
		wait
		head=$( head -1 ${outdir}/Alignment/${name}/tmpwiguniq_${name}.wig  )
		sed -i '2i'"${head}"'' ${wigFileName}
		sed -i '/	0/d' ${wigFileName}
		#generation of tdf file from wig file containing only uniquely aligned reads
		./igvtools toTDF ${wigFileName} ${tdfFileName} ${GENOME_VERSION}

	else 

		echo "Wig files generation FAILED : The file ${BIN}/${IGV_VERSION}/genomes/${GENOME_VERSION}.genome does not exist !! "

	fi
	
else 

	echo "Problem with ${outdir}/Alignment/${name}/tmpsamuniq_${name}.sam file "

fi
	
wait
#gzip tdf file
gzip ${tdfFileName}
#remove temporary files
rm ${outdir}/Alignment/${name}/tmp*
echo "`date`: Wig & tdf files generation : done "
