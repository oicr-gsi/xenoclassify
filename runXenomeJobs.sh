#!/usr/bin/env bash
for f in {1,6,11,16,21,26,31,36}
	do
	file=/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/logFiles/AMLXP_fastq_to_Align_edited.txt;
	str=$((f+1));
	str=$(cat $file | sed -n ${str}p);
	str=${str%'AMLXP'*};
	str=$(echo ${str##*'SWID'} | tr -cd '[:digit:]');
	read1=$((f+2));
	read1=$(cat $file | sed -n ${read1}p);
	read2=$((f+3));
	read2=$(cat $file | sed -n ${read2}p);
	index=/oicr/data/genomes/xenome/2013Mar/idx;

	qsub -cwd -b y -N AMLXP_Xenome_${str} -l h_vmem=30g -m beas -M "Heather.D'Souza@oicr.on.ca" "module load xenome; xenome classify -T 8 -P $index --pairs -i $read1 -i $read2 --output-filename-prefix $str";
done; 