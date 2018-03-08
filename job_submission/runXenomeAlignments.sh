#!/usr/bin/env bash

array=(
6816879_ambiguous_1.fastq
6816879_ambiguous_2.fastq
6816879_both_1.fastq
6816879_both_2.fastq
6816879_graft_1.fastq
6816879_graft_2.fastq
6816879_host_1.fastq
6816879_host_2.fastq
6816879_neither_1.fastq
6816879_neither_2.fastq );

ref=/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.12/hg19_random.fa;
animal=human;

for i in {1,2}
	do 
	if [ $i -eq 2 ]
	then 
		ref=/oicr/data/genomes/mus_musculus/UCSC/Genomic/mm10_random/bwamem/0.7.12/mm10.fa;
		animal=mouse;
	fi
	for j in {0,2,4,6,8}
		do
		name=$(echo ${array[$j]} | cut -d "_" -f 2);
		qsub -cwd -b y -N  "${animal}_${name}.bam" -m beas -M "Heather.D'Souza@oicr.on.ca" -l h_vmem=10g "module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M $ref ${array[$j]} ${array[$((j+1))]} | samtools view -Sb - >  "${animal}_${name}.bam""
	done 
done