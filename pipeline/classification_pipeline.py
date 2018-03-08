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
	human_dir = re.findall('(?<=human_dir = )([0-9a-zA-Z_/.-]*$)', variable_list[6])[0]
	mouse_dir = re.findall('(?<=mouse_dir = )([0-9a-zA-Z_/.-]*$)', variable_list[7])[0]
	index = re.findall('(?<=index = )([0-9a-zA-Z_/.-]*$)', variable_list[8])[0]
	xenome_dir = re.findall('(?<=xenome_dir = )([0-9a-zA-Z_/.-]*$)', variable_list[9])[0]
	classify_results_dir = re.findall('(?<=classify_results_dir = )([0-9a-zA-Z_/.-]*$)', variable_list[10])[0]

	return project, SWID, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, index, xenome_dir, classify_results_dir

# align reads to hg19 and mm10 index
def runAlignments(SWID, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir):
	os.system("cd {:s}; qsub -cwd -b y -N _{:s}_alignment_hg19 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}.bam; samtools sort -o {:s}_sorted.bam -n {:s}.bam\""
		.format(human_dir, SWID, human_ref, read1, read2, SWID, SWID, SWID))
	os.system("cd {:s}; qsub -cwd -b y -N _{:s}_alignment_mm10 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}.bam; samtools sort -o {:s}_sorted.bam -n {:s}.bam\""
		.format(mouse_dir, SWID, mouse_ref, read1, read2, SWID, SWID, SWID))

# run xenome on sample
def runXenome(xenome_dir, project, SWID, index, read1, read2):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_run_xenome -l h_vmem=30g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load xenome; xenome classify -T 8 -P {:s} --pairs -i {:s} -i {:s}\""
		.format(xenome_dir, SWID, SWID, index, read1, read2))

# extract, tag, concatenate, and sort read names
def extractAndTag(class_file_name, o_file, tag):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_extractAndTag_{:s} -hold_jid _{:s}_run_xenome -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load python; python ../../../scripts/extractAndTag.py {:s} {:s} {:s}\""
		.format(xenome_dir, SWID, SWID, tag, SWID, class_file_name, o_file, tag))

def concatenateAndSort(graft, host, both, ambiguous, neither, SWID):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_concatenate -hold_jid _{:s}_extractAndTag_graft,_{:s}_extractAndTag_host,_{:s}_extractAndTag_both,_{:s}_extractAndTag_ambiguous,_{:s}_extractAndTag_neither -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"cat {:s} {:s} {:s} {:s} {:s} | LC_ALL=C sort -k 1 -t ' ' > xenome_sorted.txt\""
		.format(xenome_dir, SWID, SWID, SWID, SWID, SWID, SWID, SWID, graft, host, both, ambiguous, neither))

# run classifier on sample
def runClassifier(project, SWID):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_run_classifier -hold_jid _{:s}_concatenate,_{:s}_alignment_mm10,_{:s}_alignment_hg19 -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load perl; module load samtools; perl ../../../../../scripts/classifyByAS_updated.pl {:s} {:s}\""
		.format(classify_results_dir, SWID, SWID, SWID, SWID, SWID, project, SWID))

# sort read name files
def sortClassifierOutput(SWID, class_name):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_{:s}_sorted -hold_jid _{:s}_run_classifier -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"LC_ALL=C sort -k 1 -t ' ' {:s}.txt > {:s}_sorted.txt\""
		.format(classify_results_dir, SWID, SWID, class_name, SWID, class_name, class_name))

# run verification script
def runVerification(SWID):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_verification -hold_jid _{:s}_graft_sorted,_{:s}_host_sorted,_{:s}_both_sorted,_{:s}_ambiguous_sorted,_{:s}_neither_sorted -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" \"module load python; python ../../../../../scripts/verifyClassification_updated.py ../../../../xenome/{:s}/xenome_sorted.txt graft_sorted.txt host_sorted.txt both_sorted.txt ambiguous_sorted.txt neither_sorted.txt\""
		.format(classify_results_dir, SWID, SWID, SWID, SWID, SWID, SWID, SWID, SWID))

if __name__ == '__main__':

	import os
	import sys
	import re

	checkArguments()
	project, SWID, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, index, xenome_dir, classify_results_dir = assignVariables()
	# run next two commands simultaneously
	runAlignments(SWID, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir)
	runXenome(xenome_dir, project, SWID, index, read1, read2)
	
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

