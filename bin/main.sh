# USER OPTIONS
ANALYSE_PATH="/path/to/working/directory/simu03_minidata2" #YYMMDD_verN
CONFIG_FILE="mini_config.ini"
SOFTWARE_PATH="/path/to/software/directory/new_HiFiBR2_2"

#-----------------------------------------------------------------------------------------------
# PRE-OPERATING
alias python=python3.6
cd $ANALYSE_PATH
source ${ANALYSE_PATH}/${CONFIG_FILE}

#-----------------------------------------------------------------------------------------------
# OPERATING: 
# A. DISTRIBUTING FILES: to directory named by sample 
listFile=($(find . -type f \( -name "*.fastq.gz" -o -name "*.fastq" \
	-o -name "*.fq.gz" -o -name "*.fq" \) | sort | cut -b 3-)) 	# find all fq file

	# if & elif command are the same. Diff at listFile & listFile2
if [ $sequencing_mode = "SE" ]; then
	for file in $listFile; do
		substring=$(echo $file | awk -F"_" '{print $1}')
		dirName+=(${substring})
		mkdir ${substring}
		cd ${substring}
		mkdir 1_preprocess 2_mapping 3_HiFiBR
		cd ../
		mv ${substring}* ${substring}/1_preprocess 		# cannot move directory A to directory A
		mv ${substring}/1_preprocess/$file ${substring}/1_preprocess/${substring}_all.fq
	done

elif [ $sequencing_mode = "PE" ]; then
	listFile2=$(for index in "${!listFile[@]}"; do (( index % 2 )) && echo "${listFile[index]}";done)
	for file in $listFile2; do
		substring=$(echo $file | awk -F"_" '{print $1}')
		dirName+=(${substring})
		mkdir ${substring}
		cd ${substring}
		mkdir 1_preprocess 2_mapping 3_HiFiBR
		cd ../
		mv ${substring}* ${substring}/1_preprocess 		# cannot move directory A to directory A
	done
fi


#-----------------------------------------------------------------------------------------------
# B. PROCESSING & ANALYZING
for sample in ${dirName[@]}; do
	# Prepare
	echo "-----------------------------------------------------------"
	echo "ANALYZING $sample ..."
	cd $ANALYSE_PATH/$sample/		
	cd ../

#-------------------------------------------------------------------------
# 1. 1_preprocess
	echo "* STEP 1: PREPROCESSING raw reads *"
	if [ $sequencing_mode = "PE" ]; then
		cp $SOFTWARE_PATH/script/trim_qc.sh $ANALYSE_PATH/$sample/1_preprocess/
		cd $ANALYSE_PATH/$sample/1_preprocess/

	#	Check quality before execute
		R1=$(find . -type f -name "*R1*" | cut -b 3-)
		R2=$(find . -type f -name "*R2*" | cut -b 3-)

		echo "QC. Quality control raw reads"
		fastqc --extract $R1 $R2
		mv *R1*fastqc.html *R1*fastqc.zip *R1*fastqc/
		mv *R2*fastqc.html *R2*fastqc.zip *R2*fastqc/

	#	1.1. MERGE: No need to trim-adapter (merge tool could done this)
	 	echo "** 1.1. Merged raw reads **"
		flash2 $R1 $R2 -x $maxMismatchDens -o ${sample} -M $maxOverlap --interleaved-output &> log_flash.txt

	#		Check after merge
		echo "QC. Check merged reads quality"
		fastqc --extract ${sample}.notCombined.fastq ${sample}.extendedFrags.fastq
		mv ${sample}.notCombined_fastqc.html ${sample}.notCombined_fastqc.zip ${sample}.notCombined_fastqc
		mv ${sample}.extendedFrags_fastqc.html ${sample}.extendedFrags_fastqc.zip ${sample}.extendedFrags_fastqc

	#		Add suffix "_m" to merged file
		awk '{if(NR%4==1) $0=sprintf($0"_m"); print;}' ${sample}.extendedFrags.fastq > ${sample}_merged.fq

	#	1.2. TRIM QC merged & unmerged reads
		echo "** 1.2. Get high quality merged & unmerged reads **"
		./trim_qc.sh ${sample} $trimmomatic_path
		cat ${sample}_m_highqual.fq ${sample}_nC_highqual.fq > ${sample}_all.fq

	#	1.3. Move to unnecessary
		mkdir unnecessary
		mv *fastqc/ unnecessary/

		mv ${sample}.hist ${sample}.histogram unnecessary/
		cd ../../
	fi

#-------------------------------------------------------------------------
# 2. 2_mapping (Part1)
	
	echo "* STEP 2: MAPPING high-quality reads on reference *"

#		Prepare
	cp ${ANALYSE_PATH}/${sample}/1_preprocess/${sample}_all.fq ${ANALYSE_PATH}/${sample}/2_mapping
	cp $ref_path ${ANALYSE_PATH}/${sample}/2_mapping
	ref_name=$(echo $ref_path | awk -F"/" '{print $NF}') # with extension

#		Processing
	cd $ANALYSE_PATH/$sample/2_mapping
#			Align first time - mismatch 1% (remain parameters default)
	echo "** First alignment **"
	geneious -i ${sample}_all.fq $ref_name -x $SOFTWARE_PATH/script/DSB_Junc_Align1_gap4000.optionprofile \
	-et "usePadding=False" -o ${sample}_1st.sam >> ${sample}_1stAlign_log.txt
	samtools view -H ${sample}_1st.sam > ${sample}_Aligned1.sam
	samtools view -F 4 ${sample}_1st.sam >> ${sample}_Aligned1.sam
	samtools view -H ${sample}_1st.sam > ${sample}_Unaligned1.sam
	samtools view -f 4 ${sample}_1st.sam >> ${sample}_Unaligned1.sam
#			Align second time max gap - gap 1000 (remain parameters default)	
	echo "** Second alignment **"	
	geneious -i ${sample}_Unaligned1.sam $ref_name -x $SOFTWARE_PATH/script/DSB_Junc_Align2_gap4000.optionprofile \
	-et "usePadding=False" -o ${sample}_2nd.sam >> ${sample}_2ndAlign_log.txt
	samtools view -H ${sample}_2nd.sam > ${sample}_Aligned2.sam
	samtools view -F 4 ${sample}_2nd.sam >> ${sample}_Aligned2.sam


	rm ${sample}_1st.sam ${sample}_2nd.sam
	cd ../../

#-------------------------------------------------------------------------
# 3. 2_mapping (Part 2)
	echo "* STEP 3: PRE-PROCESSING SAM file *"

#		Prepare	
	cp $SOFTWARE_PATH/script/step3* ${ANALYSE_PATH}/$sample/2_mapping

#		Processing
	cd $ANALYSE_PATH/$sample/2_mapping
	chmod +x step3_SubMain.sh
	./step3_SubMain.sh ${sample} ${ANALYSE_PATH} $cutPos1 $cutPos2 $ref_path $SOFTWARE_PATH
	cd ../../
	
#-------------------------------------------------------------------------
# 3. 2_mapping (Part 2)
 	echo "** PADDING interested SAM file **"
	cd $ANALYSE_PATH/$sample/2_mapping

#	Prepare
	cp $SOFTWARE_PATH/script/PadAnalyzeJunction2.py ${ANALYSE_PATH}/$sample/2_mapping
#	Processing
 	if [ $check_step3 == "True" ]; then
		python PadAnalyzeJunction2.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH
	elif [$check_step3 == "False" ]; then
		python PadAnalyzeJunction2.py $sam_name
	fi
	cd ../../

#-------------------------------------------------------------------------	
# 4. 3_HiFiBR
	echo "STEP 4: ANALYZING DSB Junctions"
#	Prepare
	cp $SOFTWARE_PATH/script/step4* $ANALYSE_PATH/$sample/3_HiFiBR
	cp $ANALYSE_PATH/$sample/2_mapping/${sample}_Padded.sam $ANALYSE_PATH/$sample/3_HiFiBR
	
#	Processing
	cd $ANALYSE_PATH/$sample/3_HiFiBR
	python $SOFTWARE_PATH/script/step4_alt_Hi4_5.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH $cutPos1 $cutPos2 $threads $minTimeofEvent 
	cd ../
done
	

