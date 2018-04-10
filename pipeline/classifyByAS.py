
def parse_input ():
	parser = argparse.ArgumentParser(description='Classify reads as host, graft, both, neither, or ambiguous.')
	parser.add_argument('-M', '--mouse', help='BAM file for reads aligned to mouse' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-H','--human', help='BAM file for reads aligned to human' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-O', '--output', help='Output directory for output bam and fastq files. Use "." for current working directory', type=lambda x: is_valid_directory(parser, x), required=True)
	parser.add_argument('-b', '--bam', help='Ouput BAM file with reads tagged according to assigned class', action='store_true')
	parser.add_argument('-f', '--fastq', help='Output FASTQ files', action='store_true')
	parser.add_argument('-p', '--prefix', help='Prefix for BAM and FASTQ file names', type=str)
	args = parser.parse_args()
	return parser, args 

# error handling
def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error("The file %s does not exist" % arg)
    if not re.search(r'(?<=[.])\w+',arg).group(0) == 'bam':
     	parser.error("The file %s is not a bam file" % arg)
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

def is_prefix(prefix):
	if prefix is None:
		prefix = ""
	else:
		prefix = prefix + "_"
	return prefix

def create_fastq_output(output_dir, prefix):
	# create output files

	graft_1 = open("{:s}/{:s}graft_1.fastq".format(output_dir, prefix),"w")
	host_1 = open("{:s}/{:s}host_1.fastq".format(output_dir, prefix),"w")
	both_1 = open("{:s}/{:s}both_1.fastq".format(output_dir, prefix),"w")
	ambiguous_1 = open("{:s}/{:s}ambiguous_1.fastq".format(output_dir, prefix),"w")
	neither_1 = open("{:s}/{:s}neither_1.fastq".format(output_dir, prefix),"w")

	graft_2 = open("{:s}/{:s}graft_2.fastq".format(output_dir, prefix),"w")
	host_2 = open("{:s}/{:s}host_2.fastq".format(output_dir, prefix),"w")
	both_2 = open("{:s}/{:s}both_2.fastq".format(output_dir, prefix),"w")
	ambiguous_2 = open("{:s}/{:s}ambiguous_2.fastq".format(output_dir, prefix),"w")
	neither_2 = open("{:s}/{:s}neither_2.fastq".format(output_dir, prefix),"w")

	fastq_output_1 = {
		"graft": graft_1,
		"host": host_1,
		"both": both_1,
		"ambiguous": ambiguous_1,
		"neither": neither_1,
	}

	fastq_output_2 = {
		"graft": graft_2,
		"host": host_2,
		"both": both_2,
		"ambiguous": ambiguous_2,
		"neither": neither_2,
	}

	return fastq_output_1, fastq_output_2

def create_bam_output(bam, output_dir, prefix):
	# BAMs
	out_bam = AlignmentFile("{:s}/{:s}output.bam".format(output_dir, prefix), "wb", template = bam)
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
	primary = AlignmentFile("{:s}_primary.bam".format(species), "wb", template = bam)
	secondary = AlignmentFile("{:s}_secondary.bam".format(species), "wb", template = bam)
	for read in bam.fetch(until_eof=True):
		if not read.is_secondary:
			primary.write(read)
		else:
			secondary.write(read)
	return primary, secondary

def reopen_bam(bam):
	bam.close()
	bam = AlignmentFile("human_secondary.bam", 'rb')
	return bam

def get_secondary_alignment(secondary_bam):
	secondary_alignment = next(secondary_bam)
	return secondary_alignment

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

def check_secondary_alignments(name, secondary_bam, secondary_alignment):
	secondary_alignments = []
	while secondary_alignment.query_name == name:
		secondary_alignments.append(secondary_alignment)
		try:
			secondary_alignment = next(secondary_bam)
		except StopIteration:
			break
	return secondary_alignments, secondary_alignment

def append_lists(primary_alignments, secondary_alignments):
	all_reads = primary_alignments + secondary_alignments
	return all_reads

def add_tag(read, category, output_bam):
		read.set_tag('CL',category)
		output_bam.write(read)

def convert_to_fastq(read, fastq):
	fastq.write("@{:s}\n{:s}\n+\n{:s}\n".format(read.query_name, read.query_sequence, "".join([chr(x + 33) for x in read.query_qualities])))

# calculate class percentages
def calculate_percentage(count, total_count):
	percentage = float(count)/total_count*100
	return percentage

def display_output(percentages):
	sys.stdout.write("Percentage of Reads in Each Class\n\nHost:{:.2f}\nGraft:{:.2f}\nBoth:{:.2f}\nAmbiguous:{:.2f}\nNeither:{:.2f}"
		.format(percentages[0],percentages[1],percentages[2],percentages[3],percentages[4]))

def close_file(file):
	file.close()

def remove_temp_bam(file, output_dir):
	path = output_dir + '/' + file
	try:
		os.remove(path)
	except OSError:
		print ("Error: File %s does not exist." % (path))


# output results
if __name__ == '__main__':
	
	# import required modules
	import sys
	import argparse
	import re
	import os
	import tempfile
	from pysam import AlignmentFile 

	try:
		# Python 3
		from itertools import zip_longest
	except ImportError:
		# Python 2
		from itertools import izip_longest as zip_longest

	parser, args = parse_input()
	is_prefix_dependent(parser, args)
	mouse_bam, human_bam, output_dir, is_bam, is_fastq, prefix = initialize_user_input(args)
	prefix = is_prefix(prefix)
	if is_fastq:
		fastq_output_1, fastq_output_2 = create_fastq_output(output_dir, prefix)
	if is_bam:
		output_bam = create_bam_output(human_bam, output_dir, prefix)
	counters = initialize_counters()

	# mouse_primary, mouse_secondary = filter_bam(mouse_bam, 'mouse') # change variable names to graft_bam and host_bam?
	# human_primary, human_secondary = filter_bam(human_bam, 'human')
	human_primary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/testing_classifier_scripts/test/human_primary.bam', 'rb')
	mouse_primary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/testing_classifier_scripts/test/mouse_primary.bam', 'rb')	
	human_secondary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/testing_classifier_scripts/test/human_secondary.bam', 'rb')
	mouse_secondary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/testing_classifier_scripts/test/mouse_secondary.bam', 'rb')

	human_secondary = reopen_bam(human_secondary)
	secondary_alignment = get_secondary_alignment(human_secondary)

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
			secondary_alignments, secondary_alignment = check_secondary_alignments(hr1_name, human_secondary, secondary_alignment)
			all_reads = append_lists(human_reads, secondary_alignments)
			if is_bam:
				for read in all_reads:
					add_tag(read, classification, output_bam)
			if is_fastq:
				convert_to_fastq(human_reads[0], fastq_output_1[classification])
				convert_to_fastq(human_reads[1], fastq_output_2[classification])

	# output
	for key in counters:
		percentages.append(calculate_percentage(counters[key],counters["total"]))
	display_output(percentages)

	# close files
	if is_fastq:
		for file in fastq_output_1:
			close_file(file)
		for file in fastq_output_2:
			close_file(file)
	if is_bam:
		close_file(output_bam)

	# remove temporary bams
	for file_handle, file_name in [(human_primary, 'human_primary.bam'),(human_secondary, 'human_secondary.bam'),(mouse_primary, 'mouse_primary.bam'), (mouse_secondary, 'mouse_secondary.bam')]:
		close_file(file_handle)
		remove_temp_bam(file_name, output_dir)
		