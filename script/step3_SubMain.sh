sample=$1
ANALYSE_PATH=$2
cutPos1=$3 
cutPos2=$4
ref_path=$5
SOFTWARE_PATH=$6

alias python=python3.6
# A. Fix multiple alignment info from each read
# Input: {sample}_aligned1.sam & {sample}_aligned2.sam
# Output: {sample}_fixed_aligned1.sam & {sample}_fixed_aligned2.sam
echo "** STEP 3.1: FIX SECONDARY ALIGNMENTS **"
python step3_FixAligned.py ${sample} ${ANALYSE_PATH} $SOFTWARE_PATH

# B. GET JUNCTIONS which align over 2 RE positions on reference 
# Input: {sample}_fixed_aligned1.sam & {sample}_fixed_aligned2.sam
# Output: {sample}_FilteredAligned1.sam & {sample}_FilteredAligned2.sam
echo "** STEP 3.2: GET Junctions Aligned over 2 RE positions **"
./step3_FilterMapped.sh ${sample} $cutPos1 $cutPos2


# C. ALTER & ADD alignment info
#	1. Fix NGS Error
echo "** STEP 3.3: ALTER & ADD alignment info **"
filesize=$(stat -c%s "${sample}_FilteredAligned2.sam")
echo "FILESIZE IS ${filesize}"

if [ $filesize == 0 ] # check if there is any read in 2nd mapping 
then
	mv ${sample}_FilteredAligned1.sam ${sample}_FilteredAligned.sam
	rm ${sample}_FilteredAligned2.sam
else
	samtools merge ${sample}_FilteredAligned.sam ${sample}_FilteredAligned1.sam ${sample}_FilteredAligned2.sam 
fi

python step3_FixedMappedSam2.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH
# 	2. Add End Position to sam file
python step3_Len.py ${sample}_NoNgsError.sam > ${sample}_temp.sam
samtools view -h ${sample}_temp.sam > ${sample}_len.sam
rm ${sample}_temp.sam


# D. Filter reads with full left- / right- side
python step3_Readwith2sides.py ${sample} ${ANALYSE_PATH} $ref_path $SOFTWARE_PATH $cutPos1 $cutPos2


	

