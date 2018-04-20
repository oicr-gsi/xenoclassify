# check arguments 
def checkArguments():
	if len(sys.argv) != 2:
		print("Usage: python classification_pipeline.py variables.txt")

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
	stdout = re.findall('(?<=stdout = )([0-9a-zA-Z_/.-]*$)', variable_list[11])[0]
	stderr = re.findall('(?<=stderr = )([0-9a-zA-Z_/.-]*$)', variable_list[12])[0]
	extract_script = re.findall('(?<=extract_script = )([0-9a-zA-Z_/.-]*$)', variable_list[13])[0]
	classify_script = re.findall('(?<=classify_script = )([0-9a-zA-Z_/.-]*$)', variable_list[14])[0]

	return project, SWID, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, index, xenome_dir, classify_results_dir, stdout, stderr, extract_script, classify_script

# align reads to hg19 and mm10 index
def runAlignments(project, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, stdout, stderr):
	os.system("cd {:s}; qsub -cwd -b y -N _{:s}_alignment_hg19 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}.bam; samtools sort -o {:s}_sorted.bam -n {:s}.bam\""
		.format(human_dir, project, stdout, stderr, human_ref, read1, read2, project, project, project))
	os.system("cd {:s}; qsub -cwd -b y -N _{:s}_alignment_mm10 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}.bam; samtools sort -o {:s}_sorted.bam -n {:s}.bam\""
		.format(mouse_dir, project, stdout, stderr, mouse_ref, read1, read2, project, project, project))

# run xenome on sample
def runXenome(xenome_dir, project, index, read1, read2, stdout, stderr):
	os.system("cd {:s}{:s}; pwd; qsub -cwd -b y -N _{:s}_run_xenome -l h_vmem=30g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load xenome; xenome classify -T 8 -P {:s} --pairs -i {:s} -i {:s}\""
		.format(xenome_dir, project, project, stdout, stderr, index, read1, read2))

# extract, tag, concatenate, and sort read names
def extractAndTag(class_file_name, o_file, tag, stdout, stderr, extract_script):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_extractAndTag_{:s} -hold_jid _{:s}_run_xenome -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python; python {:s} {:s} {:s} {:s}\""
		.format(xenome_dir, project, project, tag, project,  stdout, stderr, extract_script, class_file_name, o_file, tag))

def concatenateAndSort(graft, host, both, ambiguous, neither, project, stdout, stderr):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_concatenate -hold_jid _{:s}_extractAndTag_graft,_{:s}_extractAndTag_host,_{:s}_extractAndTag_both,_{:s}_extractAndTag_ambiguous,_{:s}_extractAndTag_neither -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"cat {:s} {:s} {:s} {:s} {:s} | LC_ALL=C sort -k 1 -t ' ' > xenome_sorted.txt\""
		.format(xenome_dir, project, project, project, project, project, project, project,  stdout, stderr, graft, host, both, ambiguous, neither))

# run classifier on sample
def runClassifier(project, stdout, stderr, classify_script):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_run_classifier -hold_jid _{:s}_concatenate,_{:s}_alignment_mm10,_{:s}_alignment_hg19 -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python-gsi/3.6.4; python3.6 {:s} -M {:s}{:s}_sorted.bam -H {:s}{:s}_sorted.bam -O {:s}{:s} -f\""
		.format(classify_results_dir, project, project, project, project, project,  stdout, stderr, classify_script, mouse_dir, project, human_dir, project, classify_results_dir, project))

# sort read name files
def sortClassifierOutput(project, class_name, stdout, stderr):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_{:s}_sorted -hold_jid _{:s}_run_classifier -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"LC_ALL=C sort -k 1 -t ' ' {:s}.txt > {:s}_sorted.txt\""
		.format(classify_results_dir, project, project, class_name, project,  stdout, stderr, class_name, class_name))

# run verification script
def runVerification(project, stdout, stderr):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_verification -hold_jid _{:s}_graft_sorted,_{:s}_host_sorted,_{:s}_both_sorted,_{:s}_ambiguous_sorted,_{:s}_neither_sorted -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python-gsi/3.6.4; python3.6 ../../../../../scripts/verifyClassification_updated.py ../../../../xenome/{:s}/xenome_sorted.txt graft_sorted.txt host_sorted.txt both_sorted.txt ambiguous_sorted.txt neither_sorted.txt\""
		.format(classify_results_dir, project, project, project, project, project, project, project,  stdout, stderr, project))

if __name__ == '__main__':

	import os
	import sys
	import re

	checkArguments()
	project, SWID, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, index, xenome_dir, classify_results_dir, stdout, stderr, extract_script, classify_script = assignVariables()
	# run next two commands simultaneously
	# runAlignments(project, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, stdout, stderr)
	# runXenome(xenome_dir, project, index, read1, read2, stdout, stderr)
	
	# # run next extractAndTag only after runXenome is finished
	# for tag in ['graft', 'host', 'both', 'ambiguous', 'neither']:
	# 	class_file_name = "{:s}_1.fastq".format(tag)
	# 	o_file = "xenome_{:s}.txt".format(tag)
	# 	extractAndTag(class_file_name, o_file, tag, stdout, stderr, extract_script)
	
	# # run concatenateAndSort only once all extractAndTag methods are finished
	# concatenateAndSort('xenome_graft.txt', 'xenome_host.txt', 'xenome_both.txt', 'xenome_ambiguous.txt', 'xenome_neither.txt', project, stdout, stderr)
	
	# run next method only once runAlignments() and concatenateAndSort are finished
	runClassifier(project, stdout, stderr, classify_script)
	
	# wait for runClassifier to finish
	# for class_name in ['graft_1', 'host_1', 'both_1', 'ambiguous_1', 'neither_1']:
	# 	sortClassifierOutput(project, class_name, stdout, stderr)
	
	# # wait for sortClassfierOutput to finish
	# runVerification(project, stdout, stderr)
