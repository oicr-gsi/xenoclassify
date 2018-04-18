def parse_input ():
	"""Parse command-line arguments.

	Returns:
		parser(:obj:`argparse.ArgumentParser`)
		args(:obj:`ArgumentParser.parse_args`): A namespace containing the script's input arguments.
	"""
	parser = argparse.ArgumentParser(description='Classify reads as host, graft, both, neither, or ambiguous.')
	parser.add_argument('-M', '--mouse', help='BAM file for reads aligned to mouse' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-H','--human', help='BAM file for reads aligned to human' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-O', '--output', help='Output directory for output bam and fastq files. Use "." for current working directory', type=lambda x: is_valid_directory(parser, x), required=True)
	parser.add_argument('-b', '--bam', help='Ouput BAM file with reads tagged according to assigned class', action='store_true')
	parser.add_argument('-f', '--fastq', help='Output FASTQ files', action='store_true')
	parser.add_argument('-p', '--prefix', help='Prefix for BAM and FASTQ file names', type=str)
	parser.add_argument('-n', '--neither-threshold', help='Alignment score below which all reads in a set are classified as "neither"', type=int, default=40)
	parser.add_argument('-t', '--tolerance', help='Tolerance around the mean of alignment scores for a set of reads classified as "both"', type=int, default=5)
	parser.add_argument('-d', '--difference', help='Difference between the sum of mouse and human alignment scores for a set of reads classified as "both"', type=int, default=5)
	args = parser.parse_args()
	return parser, args 

# error handling
def is_valid_file(parser, arg):
	"""Check if file is a bam and exists.
	
	Args:
		parser(:obj:`argparse.ArgumentParser`)
		arg(string): A file path.
	Returns:
		arg (string)
	Raises:
		Error if the file does not exist.
		Error if the file is not a bam file.

	"""
	if not os.path.isfile(arg):
		parser.error("The file %s does not exist" % arg)
	if not re.findall(r'(?<=[.])\w+',arg)[-1] == 'bam':
		parser.error("The file %s is not a bam file" % arg)
	else:
		return arg

def is_valid_directory(parser, arg):
	"""Check if directory exists.

	Args:
		parser(:obj:`argparse.ArgumentParser`)
		arg(string): A file path.
	Returns:
		arg (string)
	Raises:
		Error if the directory does not exist.

	"""
	if not os.path.exists(arg):
		parser.error("The output directory %s does not exist" % arg)
	else:
		if arg == '.':
			arg = os.getcwd()
		return arg

def is_prefix_dependent(parser, args):
	"""Check if 'bam' or 'fastq' are required.

	If the '--prefix' option is used, at least one of '--bam' or '--fastq' must be included as well.

	Args:
		parser(:obj:`argparse.ArgumentParser`)
		args(:obj:`ArgumentParser.parse_args`): A namespace containing the script's input arguments.
	Raises:
		Error if the '--prefix' option is used without '--bam' or '--fastq'.

	"""
	if not args.prefix is None:
		if not (args.bam is True or args.fastq is True):
			parser.error("The --prefix argument requires the use of --bam or --fastq")

def initialize_file_input(args):
	"""Set file-relate inputs to variables with concise names.

	Args:
		args (:obj:`ArgumentParser.parse_args`): A namespace containing the script's input arguments.
	Returns:
		mouse_bam(:obj:`pysam.AlignmentFile`): A bam file object of reads aligned to a mouse index.
		human_bam(:obj:`pysam.AlignmentFile`): A bam file object of reads aligned to a human index.
		output_dir(str): Directory path for output files.
		is_bam(bool): True if '--bam' specified, False if not.
		if_fastq(bool); True if '--fastq' specified, False if not.
		prefix(str or None): The string given if '--prefix' specified, 'None' if not.

	"""
	mouse_bam = args.mouse
	human_bam = args.human
	template_bam = AlignmentFile(args.human, 'rb')
	output_dir = args.output
	is_bam = args.bam
	is_fastq = args.fastq
	prefix = args.prefix
	return mouse_bam, human_bam, template_bam, output_dir, is_bam, is_fastq, prefix

def initialize_threshold_inputs(args):
	"""Set threshold-related inputs to variables with concise names.

	Args:
		args(:obj:`ArgumentParser.parse_args`): A namespace containing the script's input arguments.
	Returns:
		neither_threshold(int): Alignment score below which all reads in a set are classified as "neither."
		tolerance(int): Tolerance around the mean of alignment scores for a set of reads classified as "both."
		difference(int): Difference between the sum of mouse and human alignment scores for a set of reads classified as "both."

	"""
	neither_threshold = args.neither_threshold
	tolerance = args.tolerance
	difference = args.difference
	return neither_threshold, tolerance, difference

def is_prefix(prefix):
	"""Check if a prefix is given.

	Args:
		prefix(str or None): The string given if '--prefix' specified, 'None' if not.
	Returns:
		prefix(str): The string given prefixed by '_' if --prefix' specified, '' if not.

	"""
	if prefix is None:
		prefix = ""
	else:
		prefix = prefix + "_"
	return prefix

def create_fastq_output(output_dir, prefix):
	"""Create and wrap output fastq files for each classification.

	Args:
		output_dir(str): Directory path for output files.
		prefix(str): The string given prefixed by '_' if --prefix' specified, '' if not.
	Returns:
		fastq_output_1(dict): A wrapper for file handles to which R1 reads will be written.
		fastq_output_2(dict): A wrapper for file handles to which R2 reads will be written.

	"""
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
		"neither": neither_1
	}

	fastq_output_2 = {
		"graft": graft_2,
		"host": host_2,
		"both": both_2,
		"ambiguous": ambiguous_2,
		"neither": neither_2
	}

	return fastq_output_1, fastq_output_2

def create_bam_output(bam, output_dir, prefix):
	"""Create output bam file.

	Args:
		bam(:obj:`pysam.AlignmentFile`): A template bam file.
		output_dir(str): Directory path for output files.
		prefix(str): The string given prefixed by '_' if --prefix' specified, '' if not.
	Returns:
		out_bam (:obj:`pysam.AlignmentFile`): A bam file to which classified reads will be written.
		
	"""
	out_bam = AlignmentFile("{:s}/{:s}output.bam".format(output_dir, prefix), "wb", template = bam)
	return out_bam

def create_fastq_lists():
	"""Create and wrap empty lists for each classification.

	Returns:
		fastq_list(dict): A wrapper containing empty lists for each classification.
	"""
	graft = []
	host = []
	both = []
	ambiguous = []
	neither = []

	fastq_list = {
		"graft": graft,
		"host": host,
		"both": both,
		"ambiguous": ambiguous,
		"neither": neither
	}
	return fastq_list
	
def initialize_counters():
	"""Create and wrap counters for each classification.

	Returns:
		counters(dict): A wrapper containing counters for the number of reads belonging to each classification.

	"""
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

def filter_bam(bam_path, species):
	"""Filter primary and secondary bam alignments into separate bam files.

	Args: 
		bam(str): Path to the bam file.
		species(str): The reference genome species used to create the given bam file.
	Returns:
		primary(:obj:`pysam.AlignmentFile`): A bam file of all primary alignments.
		secondary(:obj:`pysam.AlignmentFile`): A bam file of all secondary alignments.

	"""
	bam = AlignmentFile(bam_path,"rb")
	primary = AlignmentFile("{:s}/{:s}_primary.bam".format(output_dir, species), "wb", template = bam)
	secondary = AlignmentFile("{:s}/{:s}_secondary.bam".format(output_dir, species), "wb", template = bam)
	for read in bam.fetch(until_eof=True):
		if not read.is_secondary:
			primary.write(read)
		else:
			secondary.write(read)
	bam.close()
	primary.close()
	secondary.close()

def intialize_filtered_bams(species, output_dir):
	"""Open filtered bams as read objects
	
	Args:
		species(str): The reference genome species used to create the given bam file.
		output_dir(str): Directory path for output files.
	Return:
		primary(:obj:`pysam.AlignmentFile`): A bam file containing only primary alignments.
		secondary(:obj:`pysam.AlignmentFile`): A bam file containing only secondary alignments.
	"""
	primary = AlignmentFile('{:s}/{:s}_primary.bam'.format(output_dir, species), 'rb')
	secondary = AlignmentFile('{:s}/{:s}_secondary.bam'.format(output_dir, species), 'rb')
	return primary, secondary

def get_secondary_alignment(secondary_bam):
	"""Extract an alignment from a bam file of secondary alignments.

	Args:
		secondary_bam(:obj:`pysam.AlignmentFile`): A bam file containing only secondary alignments.
	Returns:
		secondary_alignment(:obj:`AlignmentFile.AlignedSegment`): A secondary alignment.

	"""
	secondary_alignment = next(secondary_bam)
	return secondary_alignment

def get_data(read):
	"""Retrieve the name and alignment score of a read

	Args:
		read(:obj:`AlignmentFile.AlignedSegment`): An alignment.
	Returns:
		name(str): The read's query name.
		alignment_score(int): The read's alignment score.

	"""
	name = read.query_name
	alignment_score = read.get_tag('AS', with_value_type=False)
	return name, alignment_score

def check_read_names(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name):
	"""Check that all given read names match.

	Args: 
		mouse_read_1_name(str): The read name of the first read in a set of paired-reads aligned to the mouse bam.
		mouse_read_2_name(str):	The read name of the second read in a set of paried-read aligned to the mouse bam.
		human_read_1_name(str): The read name of the first read in a set of paired-reads aligned to the human bam.
		human_read_2_name(str): The read name of the second read in a set of paried-read aligned to the human bam.

	"""
	if not (mouse_read_1_name == mouse_read_2_name and human_read_1_name == human_read_2_name and mouse_read_1_name == human_read_1_name):
 		sys.exit("read names do not match: M1:{:s} M2:{:s} H1:{:s} H2:{:s}\n".format(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name))

def classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS, neither_threshold, tolerance, difference):
	"""Classify read as 'graft', 'host', 'both', 'neither', or 'ambiguous'.

	This function classifies reads based on the set of rules defined below.

	Args:
		mr1_AS(int): The alignment score of the first read in a set of paired-reads aligned to the mouse bam.
		mr2_AS(int): The alignment score of the second read in a set of paired-reads aligned to the mouse bam.
		hr1_AS(int): The alignment score of the first read in a set of paired-reads aligned to the human bam.
		hr2_AS(int): The alignment score of the second read in a set of paired-reads aligned to the human bam.
		neither_threshold(int): Alignment score below which all reads in a set are classified as "neither."
		tolerance(int): Tolerance around the mean of alignment scores for a set of reads classified as "both."
		difference(int): Difference between the sum of mouse and human alignment scores for a set of reads classified as "both."
	Returns:
		classification(str): The category to which a read has been assigned.
	"""
	scores = [mr1_AS, mr2_AS, hr1_AS, hr2_AS]
	sum_m = scores[0] + scores[1]
	sum_h = scores[2] + scores[3]
	mean = (sum_m + sum_h) / 4
	low = mean - tolerance
	high = mean + tolerance
	
	if sum_m < neither_threshold and sum_h < neither_threshold:
		classification = "neither"
	elif all(low < score < high for score in scores) or abs(sum_m - sum_h) <= difference:
		classification = "both"
	elif sum_m < sum_h:
		classification = "graft"
	elif sum_m > sum_h:
		classification = "host"
	else:
		classification = "ambiguous"

	return classification

def increment_counters(classification, counters):
	"""Increment counter associated with the current reads classification by 1.

	Args:
		classification(str): The category to which a read has been assigned.
		counters(dict): A wrapper containing counters for the number of reads belonging to each classification. 		
	"""
	counters[classification] += 1
	counters["total"] += 1

def check_secondary_alignments(name, secondary_bam, secondary_alignment):
	"""Return secondary alignments if read name matches the current primary alignment.

	Args:
		name(str): The name of the primary alignment.
		secondary(bam): A bam file containing only secondary alignments.
		secondary_alignment(:obj:`AlignmentFile.AlignedSegment`): The current secondary alignment.
	Returns:
		secondary_alignments(list): A list of secondary alignments with the same name as the primary alignment.
		secondary_alignment(:obj:`AlignmentFile.AlignedSegment`): The next alignment in the secondary bam file.

	"""
	secondary_alignments = []
	while secondary_alignment.query_name == name:
		secondary_alignments.append(secondary_alignment)
		try:
			secondary_alignment = next(secondary_bam)
		except StopIteration:
			break
	return secondary_alignments, secondary_alignment

def append_lists(primary_alignments, secondary_alignments):
	"""Concatenate the lists of primary and secondary alignments.

	Args:
		primary_alignments(list): A list of primary alignments.
		secondary_alignments(list): A list of secondary alignments.
	Returns:
		all_reads(list): A list of primary and secondary alignments.
	"""

	all_reads = primary_alignments + secondary_alignments
	return all_reads

def add_tag(read, category):
	"""Append an optional tag to an alignment.

	Args:
		read(:obj:`AlignmentFile.AlignedSegment`): An alignment.
		category(str): The category to which a read has been assigned.
	Returns:
		read(:obj:`AlignmentFile.AlignedSegment`): An alignment.
	"""

	read.set_tag('CL',category)
	return read

def convert_to_fastq(read, fastq_list):
	"""Convert the bam alignment to fastq format.

	Args:
		read(:obj:`AlignmentFile.AlignedSegment`): An alignment.
		fastq_list(list): A list of reads.

	"""
	read_fastq = "@{:s}\n{:s}\n+\n{:s}\n".format(read.query_name, read.query_sequence, "".join([chr(x + 33) for x in read.query_qualities]))
	fastq_list.append(read_fastq)

def write_to_bam(read, file):
	"""Write the pysam.AlignedSegment to a bam file.

	Args:
		read(:obj:`AlignmentFile.AlignedSegment`): An alignment.
		file(:obj:`pysam.AlignmentFile`): A bam file.

	"""
	file.write(read)

def is_list_full(list_length):
	"""Return True if the list has n elements, or False otherwise.

	Args:
		list_length(int): The length of the list.
	Returns:
		(bool) True if list is full, False otherwise.

	"""
	if list_length % 300000 == 0:
		return True
	else:
		return False

def write_to_fastq(input_list, output_file):
	"""Write all elements of a list to a fastq file.

	Args:
		input_list(list): A list of fastq reads.
		output_file(:obj:_io.TextIOWrapper): An output text file.
	
	"""
	output_file.writelines(input_list)
	del input_list[:]

def calculate_percentage(count, total_count):
	"""Calculate the percentage of reads classified in a class.

	Args:
		count(int): The number of reads in a certain category.
		total_count(int): The total number of reads across categories.
	Returns:
		percentage(float): The percentage of reads in a certain category.
	"""
	percentage = float(count)/total_count*100
	return percentage

def display_output(percentages):
	"""Display the percentage of reads classified in a class.

	Args:
		percentages(list): A list containing all percentages. 

	"""
	sys.stdout.write("Percentage of Reads in Each Class\n\nHost: {:.2f}\nGraft: {:.2f}\nBoth: {:.2f}\nAmbiguous: {:.2f}\nNeither: {:.2f}\n"
		.format(percentages[0],percentages[1],percentages[2],percentages[3],percentages[4]))

def close_file(file):
	"""Close the file.
	
	Args:
		file(:obj:pysam.AlignmentFile or :obj:_io.TextIOWrapper): A file handle.

	"""
	file.close()

def remove_temp_bam(file, output_dir):
	"""Remove the file from the file system.

	Args:
		file(:obj:pysam.AlignmentFile): A file handle.
		output_dir(str): Directory path for output files.

	"""
	path = output_dir + '/' + file
	try:
		os.remove(path)
	except OSError:
		print ("Error: File %s does not exist." % (path))

if __name__ == '__main__':
	
	# import required modules
	import sys
	import argparse
	import re
	import os
	import tempfile
	from pysam import AlignmentFile 
	import multiprocessing as mp
	import itertools

	parser, args = parse_input()
	is_prefix_dependent(parser, args)
	mouse_bam_path, human_bam_path, template_bam, output_dir, is_bam, is_fastq, prefix = initialize_file_input(args)
	neither_threshold, tolerance, difference = initialize_threshold_inputs(args)
	prefix = is_prefix(prefix)
	if is_fastq:
		fastq_output_1, fastq_output_2 = create_fastq_output(output_dir, prefix)
	if is_bam:
		output_bam = create_bam_output(template_bam, output_dir, prefix)
	counters = initialize_counters()
	
	arguments = [(mouse_bam_path, 'mouse'), (human_bam_path, 'human')]
	with mp.Pool(2) as pool:
		results = pool.starmap(filter_bam, arguments)
	mouse_primary, mouse_secondary = intialize_filtered_bams('mouse', output_dir)
	human_primary, human_secondary = intialize_filtered_bams('human', output_dir)

	# human_primary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/test_old_commit/current_changes/test_5/human_primary.bam', 'rb')
	# mouse_primary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/test_old_commit/current_changes/test_5/mouse_primary.bam', 'rb')	
	# human_secondary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/test_old_commit/current_changes/test_5/human_secondary.bam', 'rb')
	# mouse_secondary = AlignmentFile('/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/test_old_commit/current_changes/test_5/mouse_secondary.bam', 'rb')

	secondary_alignment = get_secondary_alignment(human_secondary)

	# fetch() only reads in one read at a time. read_count ensures that the following methods run only if two reads have been read in
	read_count = 0 	
	mouse_reads = [0,0]
	human_reads = [0,0]
	percentages = []
	fastq_lists_1 = create_fastq_lists()
	fastq_lists_2 = create_fastq_lists()
	# bam_list = []

	# iterate through file
	for mouse_read, human_read in itertools.zip_longest(mouse_primary.fetch(until_eof=True), human_primary.fetch(until_eof=True)):
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
			classification = classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS, neither_threshold, tolerance, difference)
			increment_counters(classification, counters)
			secondary_alignments, secondary_alignment = check_secondary_alignments(hr1_name, human_secondary, secondary_alignment)
			all_reads = append_lists(human_reads, secondary_alignments)
			if is_bam:
				for read in all_reads:
					read = add_tag(read, classification)
					write_to_bam(read, output_bam)
			if is_fastq:
				convert_to_fastq(human_reads[0], fastq_lists_1[classification])
				convert_to_fastq(human_reads[1], fastq_lists_2[classification])
				if is_list_full(counters[classification]):
					write_to_fastq(fastq_lists_1[classification], fastq_output_1[classification])
					write_to_fastq(fastq_lists_2[classification], fastq_output_2[classification])

	# write left over reads
	for current_list, current_file in itertools.zip_longest(fastq_lists_1.items(),fastq_output_1.items()):
		write_to_fastq(current_list[1], current_file[1])
	for current_list, current_file in itertools.zip_longest(fastq_lists_2.items(),fastq_output_2.items()):
		write_to_fastq(current_list[1], current_file[1])

	for key in counters:
		percentages.append(calculate_percentage(counters[key],counters["total"]))
	display_output(percentages)

	# close files
	if is_fastq:
		for key,file in fastq_output_1.items():
			close_file(file)
		for key,file in fastq_output_2.items():
			close_file(file)
	if is_bam:
		close_file(output_bam)

	for file_handle, file_name in [(human_primary, 'human_primary.bam'),(human_secondary, 'human_secondary.bam'),(mouse_primary, 'mouse_primary.bam'), (mouse_secondary, 'mouse_secondary.bam')]:
		close_file(file_handle)
		remove_temp_bam(file_name, output_dir)
		