
def parse_input ():
	parser = argparse.ArgumentParser(description='Classify reads as host, graft, both, neither, or ambiguous.')
	parser.add_argument('-M', '--mouse', help='BAM file for reads aligned to mouse' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-H','--human', help='BAM file for reads aligned to human' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-O', '--output', help='Output directory for output bam and fastq files. Use "." for current working directory', type=lambda x: is_valid_directory(parser, x), required=True)
	parser.add_argument('-b', '--bam', help='Ouput BAM file with reads tagged according to assigned class', action='store_true')
	parser.add_argument('-f', '--fastq', help='Output FASTQ files', action='store_true')
	parser.add_argument('-p', '--prefix', help='Prefix for BAM and FASTQ file names')
	args = parser.parse_args()
	print(args)
	return parser, args 

# error handling
def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("The file %s does not exist" % arg)
    # if not re.search(r'(?<=[.])\w+',arg).group(0) == 'bam':
    #  	parser.error("The file %s is not a bam file" % arg)
    else:
    	return arg

def is_valid_directory(parser, arg):
	if not os.path.exists(arg):
		parser.error("The output directory %s does not exist" % arg)
	else:
		if arg == '.':
			arg = os.getcwd()
		return arg

def is_prefix_dependent(parser, args):
	if not args.prefix is None:
		if not (args.bam is True or args.fastq is True):
			parser.error("The --prefix argument requires the use of --bam or --fastq")

def initialize_user_input(args):
	# assign arguments to script variables
	mouse_bam = AlignmentFile(args.mouse,"rb")
	human_bam = AlignmentFile(args.human,"rb")
	output_dir = args.output
	is_bam = args.bam
	is_fastq = args.fastq
	prefix = args.prefix
	return mouse_bam, human_bam, output_dir, is_bam, is_fastq, prefix

def create_fastq_output(bam):
	# create output files
	# FASTQs
	graft = AlignmentFile("graft.bam", "wb", template = bam)
	host = AlignmentFile("host.bam", "wb", template = bam)
	both = AlignmentFile("both.bam", "wb", template = bam)
	ambiguous = AlignmentFile("ambiguous.bam", "wb", template = bam)
	neither = AlignmentFile("neither.bam", "wb", template = bam)

	fastq_output = {
		"graft": graft,
		"host": host,
		"both": both,
		"ambiguous": ambiguous,
		"neither": neither,
	}
	return fastq_output

def create_bam_output(bam):
	# BAMs
	out_bam = AlignmentFile("output.bam", "wb", template = bam)
	return out_bam
	
def initialize_counters():
	# initialize counters
	host_count = 0; 
	graft_count = 0; 
	both_count = 0; 
	neither_count = 0;
	ambiguous_count = 0; 
	total_count = 0; 
	counters = {
		"host": host_count,
		"graft": graft_count,
		"both": both_count,
		"ambiguous": ambiguous_count,
		"neither": neither_count,
		"total": total_count
	}
	return counters

# filter primary and secondary alignments
def filter_bam(bam, species):
	primary = AlignmentFile("{:s}_primary.bam".format(species), "wb", template=bam)
	secondary = AlignmentFile("{:s}_secondary.bam".format(species), "wb", template=bam)
	for read in bam.fetch(until_eof=True):
		if not read.is_secondary:
			primary.write(read)
		else:
			secondary.write(read)
	return primary, secondary

def get_data(read): # check ouput of read.get_tag()
	name = read.query_name
	alignment_score = read.get_tag('AS', with_value_type=False)
	return name, alignment_score

def check_read_names(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name):
	if not (mouse_read_1_name == mouse_read_2_name and human_read_1_name == human_read_2_name and mouse_read_1_name == human_read_1_name):
 		sys.exit("read names do not match: M1:{:s} M2:{:s} H1:{:s} H2:{:s}\n".format(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name))

def classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS):
	
	scores = [mr1_AS, mr2_AS, hr1_AS, hr2_AS]
	sum_m = scores[0] + scores[1]
	sum_h = scores[2] + scores[3]
	mean = (sum_m + sum_h) / 4
	low = mean - 5
	high = mean + 5
	
	if sum_m < 40 and sum_h < 40:
		classification = "neither"
	elif all(low < score < high for score in scores) or abs(sum_m - sum_h) <= 5:
		classification = "both"
	elif sum_m < sum_h:
		classification = "graft"
	elif sum_m > sum_h:
		classification = "host"
	else:
		classification = "ambiguous"

	return classification

def increment_counters(classification, counters):
	counters[classification] += 1
	counters["total"] += 1

def add_tag(read, category, output_bam):
		read.set_tag('CL',category)
		output_bam.write(read)

def add_to_fastq(read, fastq):
	fastq.write(read)

# def activate_functions(is_bam, is_fastq):

# 	def add_tag(read, category, output_bam):
# 		read.set_tag('CL',category)
# 		output_bam.write(read)

# 	def add_to_fastq(read, fastq):
# 		fastq.write(read)

# 	if is_bam and is_fastq:
# 		def call_functions(read, category, output_bam, fastq):
# 			add_tag(read, category, output_bam)
# 			add_to_fastq(read, fastq)
# 	elif is_bam:
# 		def call_functions(read, category, output_bam):
# 			add_tag(read, category, output_bam)
# 	elif is_fastq:
# 		def call_functions(read, fastq):
# 			add_to_fastq(read, fastq)
# 	else:
# 		def call_functions():
# 			print("")

def convert_to_fastq(fastq, prefix):
	if prefix is None:
		Bedtool.bam_to_fastq(i=file, fq="{:s}_1.fastq".format(fastq), fq2="{:s}_2.fastq".format(fastq))
	else:
		Bedtool.bam_to_fastq(i=file, fq="{:s}_{:s}_1.fastq".format(prefix, fastq), fq2="{:s}_{:s}_2.fastq".format(prefix, fastq))

# calculate class percentages
def calculate_percentage(count, total_count):
	percentage = float(count)/total_count*100
	return percentage

def display_output(percentages):
	sys.stdout.write("Percentage of Reads in Each Class\n\nHost:{:.2f}\nGraft:{:.2f}\nBoth:{:.2f}\nAmbiguous:{:.2f}\nNeither:{:.2f}"
		.format(percentages[0],percentages[1],percentages[2],percentages[3],percentages[4]))

# output results
if __name__ == '__main__':
	
	# import required modules
	import sys
	import argparse
	import re
	import os
	import tempfile
	from pysam import AlignmentFile 
	from pybedtools import Bedtool
	
	try:
		# Python 3
		from itertools import zip_longest
	except ImportError:
		# Python 2
		from itertools import izip_longest as zip_longest

	parser, args = parse_input()
	is_prefix_dependent(parser, args)
	mouse_bam, human_bam, output_dir, is_bam, is_fastq, prefix = initialize_user_input(args)
	# activate_functions(is_bam, is_fastq)
	if is_fastq:
		fastq_output = create_fastq_output(human_bam)
	if is_bam:
		output_bam = create_bam_output(human_bam)
	counters = initialize_counters()

	mouse_primary, mouse_secondary = filter_bam(mouse_bam, 'mouse') # change variable names to graft_bam and host_bam?
	human_primary, human_secondary = filter_bam(human_bam, 'human')
	# human_primary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/classifyByAS_results/new_rules/6816887/hg19_primary.bam', 'rb')
	# mouse_primary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/classifyByAS_results/new_rules/6816887/mm10_primary.bam', 'rb')	
	
	# fetch() only reads in one read at a time. read_count ensures that the following methods run only if two reads have been read in
	read_count = 0 	
	mouse_reads = [0,0]
	human_reads = [0,0]
	percentages = []

	# iterate through file
	for mouse_read, human_read in zip_longest(mouse_primary.fetch(until_eof=True), human_primary.fetch(until_eof=True)):
		mouse_reads[read_count] = mouse_read
		human_reads[read_count] = human_read
		read_count += 1
		if read_count % 2 == 0:
			read_count = 0
			mr1_name, mr1_AS = get_data(mouse_reads[0])
			mr2_name, mr2_AS = get_data(mouse_reads[1])
			hr1_name, hr1_AS = get_data(human_reads[0])
			hr2_name, hr2_AS = get_data(human_reads[1])
			check_read_names(mr1_name, mr2_name, hr1_name, hr2_name)
			classification = classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS)
			increment_counters(classification, counters)
			for read in human_reads:
				if is_bam:
					add_tag(read, classification, output_bam)
				if is_fastq:
					add_to_fastq(read, fastq_output[classification])
				# call_functions(read, classification, output_bam, read, fastq_output[classification])

	# output
	for key in counters:
		percentages.append(calculate_percentage(counters[key],counters["total"]))
	display_output(percentages)

	#fastq
	if is_fastq is True:
		for file in fastq_output:
			convert_to_fastq(file, prefix)
