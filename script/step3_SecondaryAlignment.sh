#!/bin/bash
sample=$1
i=$2
wdir=$3
exp=${sample: -1}

cd ${wdir}
samtools view ${sample}_Aligned${i}.sam | awk -F"\t" '{print $1}'| sort | uniq -u > uniqNames${exp}.txt
samtools view ${sample}_Aligned${i}.sam | awk -F"\t" '{print $1}'| sort | uniq -d > secondaryAlignmentNames${exp}.txt

#samtools view -c ${sample}_a_align${i}.sam
wc -l uniqNames${exp}.txt
wc -l secondaryAlignmentNames${exp}.txt

samtools view -H ${sample}_Aligned${i}.sam > ${sample}_fixed_a_align${i}.sam
samtools view ${sample}_Aligned${i}.sam | rg -f uniqNames${exp}.txt >> ${sample}_fixed_a_align${i}.sam

samtools view -H ${sample}_Aligned${i}.sam > ${sample}_secondaryAlignment${i}.sam
samtools view ${sample}_Aligned${i}.sam | rg -f secondaryAlignmentNames${exp}.txt >> ${sample}_secondaryAlignment${i}.sam

rm uniqNames${exp}.txt secondaryAlignmentNames${exp}.txt




