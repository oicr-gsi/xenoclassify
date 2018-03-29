# check number of arguments passed is correct
# def checkArguments():
# 	if len(sys.argv) != 4:
# 		print("Usage: python classifyByAS.py mouse.bam human.bam ouput_dir")
# 		sys.exit()

def parse_input ():
	parser = argparse.ArgumentParser(description='Classify reads as host, graft, both, neither, or ambiguous.')
	parser.add_argument('-M', '--mouse', help='Bam file for reads aligned to mouse' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-H','--human', help='Bam file for reads aligned to human' , type=lambda x: is_valid_file(parser, x), required=True)
	parser.add_argument('-O', '--output', help='Output directory for output bam and fastq files. Use "." for current working directory', type=lambda x: is_valid_directory(parser, x), required=True)
	args = parser.parse_args()
	print(args)
	return args 

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

def initialize_variables(args):
	# assign arguments to script variables
	mouse_bam = args.mouse
	human_bam = args.human
	output_dir = args.output

	# create output files
	graft = open("{:s}/graft.txt".format(output_dir),"w")
	host = open("{:s}/host.txt".format(output_dir),"w")
	both = open("{:s}/both.txt".format(output_dir),"w")
	ambiguous = open("{:s}/ambiguous.txt".format(output_dir),"w")
	neither = open("{:s}/neither.txt".format(output_dir),"w")

	# initialize counters
	host_count = 0; 
	graft_count = 0; 
	both_count = 0; 
	neither_count = 0;
	ambiguous_count = 0; 
	total_count = 0; 
	return (mouse_bam, human_bam, output_dir, graft, host, both, ambiguous, neither, host_count, 
	 	graft_count, both_count, ambiguous_count, neither_count, total_count)

# filter primary and secondary alignments
def filter_bam(bam, species):
	sam = pysam.AlignmentFile(bam, "rb")
	primary = pysam.AlignmentFile("{:s}_primary.bam".format(species), "wb", template=sam)
	secondary = pysam.AlignmentFile("{:s}_secondary.bam".format(species), "wb", template=sam)
	count_primary=0
	count_secondary=0	
	for read in sam.fetch(until_eof=True):
		if not read.is_secondary:
			primary.write(read)
			count_primary += 1
		else:
			secondary.write(read)
			count_secondary += 1
	return sam, primary, secondary, count_primary, count_secondary

def get_data(read): # check ouput of read.get_tag()
	name = read.query_name
	alignment_score = read.get_tag('AS', with_value_type=False)
	cigar = read.cigarstring
	return name, alignment_score, cigar

def check_read_names(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name):
	if not (mouse_read_1_name == mouse_read_2_name and human_read_1_name == human_read_2_name and mouse_read_1_name == human_read_1_name):
 		sys.exit("read names do not match: M1:{:s} M2:{:s} H1:{:s} H2:{:s}\n".format(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name))

def classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS, name, host, host_count, graft, graft_count,
	both, both_count, ambiguous, ambiguous_count, neither, neither_count,  total_count, m1_cigar, m2_cigar, hr1_cigar, hr2_cigar):
	
	scores = [mr1_AS, mr2_AS, hr1_AS, hr2_AS]
	sum_m = scores[0] + scores[1]
	sum_h = scores[2] + scores[3]
	mean = (sum_m + sum_h) / 4
	low = mean - 5
	high = mean + 5
	
	if sum_m < 40 and sum_h < 40:
		neither.write('{:s}\n'.format(name))
		neither_count += 1
	elif all(low < score < high for score in scores) or abs(sum_m - sum_h) <= 5:
		both.write('{:s}\n'.format(name))
		both_count += 1
	elif sum_m < sum_h:
		graft.write('{:s}\n'.format(name))
		graft_count += 1
	elif sum_m > sum_h:
		host.write('{:s}\n'.format(name))
		host_count += 1
	else:
		ambiguous.write('{:s}\n'.format(name))
		ambiguous_count += 1
	total_count += 1
	return graft_count, host_count, both_count, ambiguous_count, neither_count, total_count

# calculate class percentages
def calculate_percentages(graft_count, host_count, both_count, ambiguous_count, neither_count, total_count):
	percentages = ['','','','','']
	percentages[0] = float(graft_count)/total_count*100
	percentages[1] = float(host_count)/total_count*100
	percentages[2] = float(both_count)/total_count*100
	percentages[3] = float(ambiguous_count)/total_count*100
	percentages[4] = float(neither_count)/total_count*100
	return percentages

def display_output(percentages):
	sys.stdout.write("Percentage of Reads in Each Class\n\nGraft:{:.2f}\nHost:{:.2f}\nBoth:{:.2f}\nAmbiguous:{:.2f}Neither:{:.2f}\n"
		.format(percentages[0],percentages[1],percentages[2],percentages[3],percentages[4]))

# output results
if __name__ == '__main__':
	
	# import required modules
	import pysam 
	import sys
	import argparse
	import re
	import os
	
	try:
		# Python 3
		from itertools import zip_longest
	except ImportError:
		# Python 2
		from itertools import izip_longest as zip_longest

	# checkArguments()
	args = parse_input()
	(mouse_bam, human_bam, output_dir, graft, host, both, ambiguous, neither, host_count, 
	 	graft_count, both_count, ambiguous_count, neither_count, total_count) = initialize_variables(args)
	mouse_sam, mouse_primary, mouse_secondary, count_primary, count_secondary = filter_bam(mouse_bam, 'mouse') # change variable names to graft_bam and host_bam?
	human_sam, human_primary, human_secondary, count_primary, count_secondary = filter_bam(human_bam, 'human')
	
	# fetch() only reads in one read at a time. read_count ensures that the following methods run only if two reads have been read in
	read_count = 0 	
	reads_mouse = [0,0]
	reads_human = [0,0]

	# iterate through file
	for read_mouse, read_human in zip_longest(mouse_primary.fetch(until_eof=True), human_primary.fetch(until_eof=True)):
		reads_mouse[read_count] = read_mouse
		reads_human[read_count] = read_human
		read_count += 1
		if read_count % 2 == 0:
			read_count = 0
			mr1_name, mr1_AS, mr1_cigar = get_data(reads_mouse[0])
			mr2_name, mr2_AS, mr2_cigar = get_data(reads_mouse[1])
			hr1_name, hr1_AS, hr1_cigar = get_data(reads_human[0])
			hr2_name, hr2_AS, hr2_cigar = get_data(reads_human[1])
			check_read_names(mr1_name, mr2_name, hr1_name, hr2_name)
			graft_count, host_count, both_count, ambiguous_count, neither_count, total_count = classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS, 
 			hr1_name, host, host_count, graft, graft_count, both, both_count, ambiguous, ambiguous_count, neither, neither_count, total_count, mr1_cigar, mr2_cigar, hr1_cigar, hr2_cigar)

	percentages = calculate_percentages(graft_count, host_count, both_count, ambiguous_count, neither_count, total_count)
	display_output(percentages)
