# check arguments 
def checkArguments():
	if len(sys.argv) != 2:
		print("Usage: python classification_pipeline.py variables.txt")

def assignVariables():
	variable_list = [line.rstrip('\n') for line in open(sys.argv[1])]
	sample = re.findall('(?<=sample = )([0-9a-zA-Z_/.-]*$)', variable_list[0])[0]
	human_ref = re.findall('(?<=human_ref = )([0-9a-zA-Z_/.-]*$)', variable_list[1])[0]
	mouse_ref = re.findall('(?<=mouse_ref = )([0-9a-zA-Z_/.-]*$)', variable_list[2])[0]
	read1 = re.findall('(?<=read1 = )([0-9a-zA-Z_/.-]*$)', variable_list[3])[0]
	read2 = re.findall('(?<=read2 = )([0-9a-zA-Z_/.-]*$)', variable_list[4])[0]
	human_bam = re.findall('(?<=human_bam = )([0-9a-zA-Z_/.-]*$)', variable_list[5])[0]
	mouse_bam = re.findall('(?<=mouse_bam = )([0-9a-zA-Z_/.-]*$)', variable_list[6])[0]
	index = re.findall('(?<=index = )([0-9a-zA-Z_/.-]*$)', variable_list[7])[0]
	xenome_dir = re.findall('(?<=xenome_dir = )([0-9a-zA-Z_/.-]*$)', variable_list[8])[0]
	classify_results_dir = re.findall('(?<=classify_results_dir = )([0-9a-zA-Z_/.-]*$)', variable_list[9])[0]
	stdout = re.findall('(?<=stdout = )([0-9a-zA-Z_/.-]*$)', variable_list[10])[0]
	stderr = re.findall('(?<=stderr = )([0-9a-zA-Z_/.-]*$)', variable_list[11])[0]
	extract_script = re.findall('(?<=extract_script = )([0-9a-zA-Z_/.-]*$)', variable_list[12])[0]
	classify_script = re.findall('(?<=classify_script = )([0-9a-zA-Z_/.-]*$)', variable_list[13])[0]
	verify_script = re.findall('(?<=verify_script = )([0-9a-zA-Z_/.-]*$)', variable_list[14])[0]
	extract_names_script = re.findall('(?<=extract_names_script = )([0-9a-zA-Z_/.-]*$)', variable_list[15])[0]

	return sample, human_ref, mouse_ref, read1, read2, human_bam, mouse_bam, index, xenome_dir, classify_results_dir, stdout, stderr, extract_script, classify_script, verify_script, extract_names_script

# align reads to hg19 and mm10 index
def runAlignments(sample, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, stdout, stderr):
	os.system("cd {:s}; qsub -cwd -b y -N _{:s}_alignment_hg19 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}.bam; samtools sort -o {:s}_sorted.bam -n {:s}.bam\""
		.format(human_dir, sample, stdout, stderr, human_ref, read1, read2, sample, sample, sample))
	os.system("cd {:s}; qsub -cwd -b y -N _{:s}_alignment_mm10 -l h_vmem=15g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load bwa/0.7.12; module load samtools; bwa mem -t 8 -M {:s} {:s} {:s} | samtools view -Sb - > {:s}.bam; samtools sort -o {:s}_sorted.bam -n {:s}.bam\""
		.format(mouse_dir, sample, stdout, stderr, mouse_ref, read1, read2, sample, sample, sample))

# run xenome on sample
def runXenome(xenome_dir,sample, index, read1, read2, stdout, stderr):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_run_xenome -l h_vmem=30g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load xenome; xenome classify -T 8 -P {:s} --pairs -i {:s} -i {:s}\""
		.format(xenome_dir,sample,sample, stdout, stderr, index, read1, read2))

# extract, tag, concatenate, and sort read names
def extractAndTag(sample, class_file_name, o_file, tag, stdout, stderr, extract_script):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_extractAndTag_{:s} -hold_jid _{:s}_run_xenome -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python; python {:s} {:s} {:s} {:s}\""
		.format(xenome_dir, sample, sample, tag, sample,  stdout, stderr, extract_script, class_file_name, o_file, tag))

def concatenateAndSort(graft, host, both, ambiguous, neither, sample, stdout, stderr):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_concatenate -hold_jid _{:s}_extractAndTag_graft,_{:s}_extractAndTag_host,_{:s}_extractAndTag_both,_{:s}_extractAndTag_ambiguous,_{:s}_extractAndTag_neither -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"cat {:s} {:s} {:s} {:s} {:s} | LC_ALL=C sort -k 1 -t ' ' > xenome_sorted.txt\""
		.format(xenome_dir, sample, sample, sample, sample, sample, sample, sample,  stdout, stderr, graft, host, both, ambiguous, neither))

# run classifier on sample
def runClassifier(classify_results_dir, sample, stdout, stderr, classify_script, human_bam, mouse_bam):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_run_classifier -hold_jid _{:s}_concatenate -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python-gsi/3.6.4; python3.6 {:s} -H {:s} -G {:s} -O {:s}{:s} -f -b\""
		.format(classify_results_dir, sample, sample, sample, stdout, stderr, classify_script, mouse_bam, human_bam, classify_results_dir, sample))

def extractNames(classify_results_dir, sample, extract_names_script, in_file, out_file, tag):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_extract_read_names_{:s} -hold_jid _{:s}_run_classifier -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python-gsi/3.6.4; python3.6 {:s} {:s} {:s}\""
		.format(classify_results_dir, sample, sample, tag, sample, stdout, stderr, extract_names_script, in_file, out_file))
	
# sort read name files
def sortClassifierOutput(classify_results_dir, sample, class_name, stdout, stderr):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_{:s}_sorted -hold_jid _{:s}_extract_read_names_{:s} -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"LC_ALL=C sort -k 1 -t ' ' {:s}_1.txt > {:s}_sorted.txt\""
		.format(classify_results_dir, sample, sample, class_name, sample, class_name,  stdout, stderr, class_name, class_name))

# run verification script
def runVerification(sample, stdout, stderr, verify_script):
	os.system("cd {:s}{:s}; qsub -cwd -b y -N _{:s}_verification -hold_jid _{:s}_graft_sorted,_{:s}_host_sorted,_{:s}_both_sorted,_{:s}_neither_sorted -l h_vmem=10g -m beas -M \"Heather.D'Souza@oicr.on.ca\" -o {:s} -e {:s} \"module load python-gsi; python {:s} {:s}{:s}/xenome_sorted.txt graft_sorted.txt host_sorted.txt both_sorted.txt neither_sorted.txt\""
		.format(classify_results_dir, sample, sample,sample,sample,sample,sample, stdout, stderr, verify_script, xenome_dir, sample))

if __name__ == '__main__':

	import os
	import sys
	import re

	checkArguments()
	sample, human_ref, mouse_ref, read1, read2, human_bam, mouse_bam, index, xenome_dir, classify_results_dir, stdout, stderr, extract_script, classify_script, verify_script, extract_names_script = assignVariables()
	# run next two commands simultaneously
	# runAlignments(project, human_ref, mouse_ref, read1, read2, human_dir, mouse_dir, stdout, stderr)
	runXenome(xenome_dir, sample, index, read1, read2, stdout, stderr)
	
	# # # run next extractAndTag only after runXenome is finished
	for tag in ['graft', 'host', 'both', 'ambiguous', 'neither']:
		class_file_name = "{:s}_1.fastq".format(tag)
		o_file = "xenome_{:s}.txt".format(tag)
		extractAndTag(sample, class_file_name, o_file, tag, stdout, stderr, extract_script)
	
	# # # run concatenateAndSort only once all extractAndTag methods are finished
	concatenateAndSort('xenome_graft.txt', 'xenome_host.txt', 'xenome_both.txt', 'xenome_ambiguous.txt', 'xenome_neither.txt', sample, stdout, stderr)
	
	# run next method only once runAlignments() and concatenateAndSort are finished
	runClassifier(classify_results_dir, sample, stdout, stderr, classify_script, human_bam, mouse_bam)

	for in_file, out_file, tag in [('graft_1.fastq', 'graft_1.txt','graft'),('host_1.fastq', 'host_1.txt','host'),('both_1.fastq', 'both_1.txt','both'),('neither_1.fastq', 'neither_1.txt','neither')]:
		extractNames(classify_results_dir, sample, extract_names_script, in_file, out_file, tag)
	
	# # wait for runClassifier to finish
	for class_name in ['graft', 'host', 'both', 'neither']:
		sortClassifierOutput(classify_results_dir, sample, class_name, stdout, stderr)
	
	# # wait for sortClassfierOutput to finish
	runVerification(sample, stdout, stderr, verify_script)
