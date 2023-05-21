sample=$1
cutPos1=$2
cutPos2=$3

for i in 1 2
do
	# Add YS tag: right-most length for read mapped onto ref
	python3 step3_Len.py ${sample}_fixed_a_align${i}.sam > alt_Contig${i}.sam

	# sam to bam file
	samtools view -bS alt_Contig${i}.sam > alt_Contig${i}.bam

	# Filter right-most pos if > first cut position
	bamtools filter -in alt_Contig${i}.bam -out filter_right${i}.bam -tag "YS:>${cutPos2}"

	# Filter left-most pos if < second cut position
	samtools view -HS filter_right${i}.bam > ${sample}_FilteredAligned${i}.sam
	samtools view -S filter_right${i}.bam | awk -F"\t" -v var="${cutPos1}" '$4<=var {print $0}' >> ${sample}_FilteredAligned${i}.sam 

	rm alt_Contig${i}.bam filter_right${i}.bam alt_Contig${i}.sam
done
