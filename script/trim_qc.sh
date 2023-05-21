# VARIABLES
sample=$1
trimmomatic_path=$2

java -jar $trimmomatic_path SE ${sample}_merged.fq ${sample}_m_highqual.fq SLIDINGWINDOW:4:15 MINLEN:150 2> log_m_qc.txt

java -jar $trimmomatic_path SE ${sample}.notCombined.fastq ${sample}_nC_highqual.fq SLIDINGWINDOW:4:15 MINLEN:150 2> log_nC_qc.txt

fastqc --extract ${sample}_m_highqual.fq ${sample}_nC_highqual.fq
mv ${sample}_m_highqual_fastqc.html ${sample}_m_highqual_fastqc.zip ${sample}_m_highqual_fastqc
mv ${sample}_nC_highqual_fastqc.html ${sample}_nC_highqual_fastqc.zip ${sample}_nC_highqual_fastqc
