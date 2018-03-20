# check number of arguments passed is correct
def checkArguments():
	if len(sys.argv) != 4:
		print("Usage: python classifyByAS.py mouse.bam human.bam ouput_dir")
		sys.exit()

def initializeVariables():
	# assign arguments to script variables
	mouse_bam = sys.argv[1]
	human_bam = sys.argv[2]
	output_dir = sys.argv[3]

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

	return (mouse_bam, human_bam, output_dir, graft, host, both, ambiguous, neither,
		host_count, graft_count, both_count, neither_count, ambiguous_count, total_count)

# filter primary and secondary alignments
def filter_bam(bam):
	cmd_header = ['samtools', 'view', '-H', bam]
	cmd_primary = ['samtools', 'view', '-F', '256', bam]
	cmd_secondary = ['samtools', 'view', '-f', '256', bam]
	header = subprocess.Popen(cmd_header, stdout=subprocess.PIPE).stdout.readline
	primary = subprocess.Popen(cmd_primary, stdout=subprocess.PIPE).stdout.readline 
	secondary = subprocess.Popen(cmd_secondary, stdout=subprocess.PIPE).stdout.readline
	return header, primary, secondary

def getData(read): # check ouput of read.get_tag()
	read = re.sub('\s+',' ',read).strip()
	name = read.split(' ')[0].strip()
	alignment_score = int(re.findall('(?<=AS:i:)([0-9]*)', read)[0])
	return name, alignment_score

def checkReadNames(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name):
	if not (mouse_read_1_name == mouse_read_2_name and human_read_1_name == human_read_2_name and mouse_read_1_name == human_read_1_name):
 		sys.exit("read names do not match: M1:{:s} M2:{:s} H1:{:s} H2:{:s}\n".format(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name))

def classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS, name, host, host_count, graft, graft_count,
	both, both_count, neither, neither_count, ambiguous, ambiguous_count, total_count):
	if mr1_AS < 104 and mr2_AS < 104 and (hr1_AS > 104 or hr2_AS > 104):
		graft.write('{:s}\n'.format(name))
		graft_count += 1
	elif (mr1_AS > 104 or mr2_AS > 104) and hr1_AS < 104 and hr2_AS < 104:
		host.write('{:s}\n'.format(name))
		host_count += 1
	elif mr1_AS > 104 and mr2_AS > 104 and hr1_AS > 104 and hr2_AS > 104:
		both.write('{:s}\n'.format(name))
		both_count += 1
	elif mr1_AS < 10 and mr2_AS < 10 and hr1_AS < 10 and hr2_AS < 10:
		neither.write('{:s}\n'.format(name))
		neither_count += 1
	else:
		ambiguous.write('{:s}\n'.format(name))
		ambiguous_count += 1
	total_count += 1
	return graft_count, host_count, both_count, neither_count, ambiguous_count, total_count

# calculate class percentages
def calculatePercentages(graft_count, host_count, both_count, ambiguous_count, neither_count, total_count):
	percentages = ['','','','','']
	percentages[0] = float(graft_count)/total_count
	percentages[1] = float(host_count)/total_count
	percentages[2] = float(both_count)/total_count
	percentages[3] = float(ambiguous_count)/total_count
	percentages[4] = float(neither_count)/total_count
	return percentages

def displayOutput(percentages):
	sys.stdout.write("Percentage of Reads in Each Class\n\nGraft:{:s}\nHost:{:s}\nBoth:{:s}\nNeither:{:s}\nAmbiguous:{:s}"
		.format(percentages[0],percentages[1],percentages[2],percentages[3],percentages[4],))

# output results
if __name__ == '__main__':
	
	# import required modules
	import subprocess
	import sys
	import argparse
	import re
	
	try:
    	# Python 3
		from itertools import zip_longest
	except ImportError:
    	# Python 2
		from itertools import izip_longest as zip_longest

	checkArguments()
	(mouse_bam, human_bam, output_dir, graft, host, both, ambiguous, neither, host_count, 
		graft_count, both_count, neither_count, ambiguous_count, total_count) = initializeVariables()
	
	mouse_header, mouse_primary, mouse_secondary = filter_bam(mouse_bam) # change variable names to graft_bam and host_bam?
	human_header, human_primary, human_secondary = filter_bam(human_bam)
	
	read_count = 0 	
	reads_mouse = [0,0]
	reads_human = [0,0]
	
	# iterate through file
	for read_mouse, read_human in zip_longest(iter(mouse_primary,''), iter(human_primary,'')):
		reads_mouse[read_count] = read_mouse.decode('utf-8')
		reads_human[read_count] = read_human.decode('utf-8')
		read_count += 1
		if read_count % 2 == 0:
			read_count = 0
			mr1_name, mr1_AS = getData(reads_mouse[0])
			mr2_name, mr2_AS = getData(reads_mouse[1])
			hr1_name, hr1_AS = getData(reads_human[0])
			hr2_name, hr2_AS = getData(reads_human[1])
			checkReadNames(mr1_name, mr2_name, hr1_name, hr2_name)
			graft_count, host_count, both_count, neither_count, ambiguous_count, total_count = classify(mr1_AS, mr2_AS, hr1_AS, hr2_AS, 
 	        hr1_name, host, host_count, graft, graft_count, both, both_count, neither, neither_count, ambiguous, ambiguous_count, total_count)

	percentages = calculatePercentages(graft_count, host_count, both_count, ambiguous_count, neither_count, total_count)
	displayOutput(percentages)
	