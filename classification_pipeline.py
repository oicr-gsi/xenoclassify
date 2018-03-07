# check arguments 
def checkArguments():
	if len(sys.argv) != 2:
		print "Usage: python classification_pipeline.py variables.txt"

def assignVariables():
	variable_list = [line.rstrip('\n') for line in open(sys.argv[1])]
	project = re.findall('(?<=project = )([0-9a-zA-Z_/]*$)', variable_list[0])[0]
	SWID = re.findall('(?<=SWID = )([0-9a-zA-Z_/]*$)', variable_list[1])[0]
	human_ref = re.findall('(?<=human_ref = )([0-9a-zA-Z_/.-]*$)', variable_list[2])[0]
	mouse_ref = re.findall('(?<=mouse_ref = )([0-9a-zA-Z_/.-]*$)', variable_list[3])[0]
	read1 = re.findall('(?<=read1 = )([0-9a-zA-Z_/.-]*$)', variable_list[4])[0]
	read2 = re.findall('(?<=read2 = )([0-9a-zA-Z_/.-]*$)', variable_list[5])[0]
	o_file_human = re.findall('(?<=o_file_human = )([0-9a-zA-Z_/.-]*$)', variable_list[6])[0]
	o_file_mouse = re.findall('(?<=o_file_mouse = )([0-9a-zA-Z_/.-]*$)', variable_list[7])[0]
	index = re.findall('(?<=index = )([0-9a-zA-Z_/.-]*$)', variable_list[8])[0]
	xenome_directory = re.findall('(?<=xenome_directory = )([0-9a-zA-Z_/.-]*$)', variable_list[9])[0]

	return project, SWID, human_ref, mouse_ref, read1, read2, o_file_human, o_file_mouse, index, xenome_directory

# align reads to hg19 and mm10 index
def runAlignments(SWID, human_ref, mouse_ref, read1, read2, o_file_human, o_file_mouse):
	os.system("qsub -cwd -b y -N _alignment_hg19 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}; samtools sort -o {:s}_sorted.bam -n {:s}\""
		.format(human_ref, read1, read2, o_file_human, SWID, o_file_human))
	os.system("qsub -cwd -b y -N _alignment_mm10 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}; samtools sort -o {:s}_sorted.bam -n {:s}\""
		.format(mouse_ref, read1, read2, o_file_mouse, SWID, o_file_mouse))

# run xenome on sample
def runXenome(xenome_directory, project, SWID, index, read1, read2):
	os.system("cd {:s}; qsub -cwd -b y -N _run_xenome -l h_vmem=30g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load xenome; xenome classify -T 8 -P {:s} --pairs -i {:s} -i {:s}\""
		.format(xenome_directory, index, read1, read2))

# extract, tag, concatenate, and sort read names
def extractAndTag(class_file_name, o_file, tag):
	os.system("qsub -cwd -b y -N _extractAndTag_{:s} -hold_jid _run_xenome -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load python; python extractAndTag.py {:s} {:s} {:s}\""
		.format(tag, class_file_name, o_file, tag))

def concatenateAndSort(graft, host, both, ambiguous, neither, SWID):
	os.system("qsub -cwd -b y -N _concatenate -hold_jid _extractAndTag_graft,_extractAndTag_host,_extractAndTag_both,_extractAndTag_ambiguous,_extractAndTag_neither -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"cat {:s} {:s} {:s} {:s} {:s} | LC_ALL=C sort -k 1 -t ' ' > {:s}_xenome_sorted.txt\""
		.format(graft, host, both, ambiguous, neither, SWID))

# run classifier on sample
def runClassifier(project, SWID):
	os.system("qsub -cwd -b y -N _run_classifier -hold_jid _concatenate,_alignment_mm10,_alignment_hg19 -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load python; python classifyByAS_updated.py {:s} {:s}\""
		.format(project, SWID))

# sort read name files
def sortClassifierOutput(SWID, class_name):
	os.system("qsub -cwd -b y -N _{:s}_sorted -hold_jid _run_classifier -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"LC_ALL=C sort -k 1 -t ' ' {:s}.txt > {:s}_{:s}_sorted.txt\""
		.format(class_name, class_name, SWID, class_name))

# run verification script
def runVerification(SWID):
	os.system("qsub -cwd -b y -N _{:s}_verification -hold_jid _graft_sorted,_host_sorted,_both_sorted,_ambiguous_sorted,_neither_sorted -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load python; python verifyClassification_updated.py {:s}_xenome_sorted.txt {:s}_graft_sorted.txt {:s}_host_sorted.txt {:s}_both_sorted.txt {:s}_ambiguous_sorted.txt {:s}_neither_sorted.txt\""
		.format(SWID, SWID, SWID, SWID, SWID, SWID, SWID))

if __name__ == '__main__':

	import os
	import sys
	import re

	checkArguments()
	project, SWID, human_ref, mouse_ref, read1, read2, o_file_human, o_file_mouse, index, xenome_directory = assignVariables()
	# run next two commands simultaneously
	runAlignments(SWID, human_ref, mouse_ref, read1, read2, o_file_human, o_file_mouse)
	runXenome(xenome_directory, project, SWID, index, read1, read2)
	
	# run next extractAndTag only after runXenome is finished
	for tag in ['graft', 'host', 'both', 'ambiguous', 'neither']:
		class_file_name = "{:s}_1.fastq".format(tag)
		o_file = "xenome_{:s}.txt".format(tag)
		extractAndTag(class_file_name, o_file, tag)
	
	# run concatenateAndSort only once all extractAndTag methods are finished
	concatenateAndSort('xenome_graft.txt', 'xenome_host.txt', 'xenome_both.txt', 'xenome_ambiguous.txt', 'xenome_neither.txt', SWID)
	
	# run next method only once runAlignments() and concatenateAndSort are finished
	runClassifier(project, SWID)
	
	# wait for runClassifier to finish
	for class_name in ['graft', 'host', 'both', 'ambiguous', 'neither']:
		sortClassifierOutput(SWID, class_name)
	
	# wait for sortClassfierOutput to finish
	runVerification(SWID)
