# check number of arguments passed is correct
def checkArguments():
	if len(sys.argv) != 2:
		print "Usage: python classifyByAS.py mouse.bam human.bam ouput_dir"
		sys.exit()

def initializeVariables():
	# assign arguments to script variables
	mouse_bam=sys.argv[1]
	human_bam=sys.argv[2]
	# create output files
	graft = open("{:s}/graft.txt".format(sys.argv[3]),"w")
	host = open("{:s}/host.txt".format(sys.argv[3]),"w")
	both = open("{:s}/both.txt".format(sys.argv[3]),"w")
	ambiguous = open("{:s}/ambiguous.txt".format(sys.argv[3]),"w")
	neither = open("{:s}/neither.txt".format(sys.argv[3]),"w")
	# initialize counters
	total_count = 0; 
	host_count = 0; 
	graft_count = 0; 
	both_count = 0; 
	neither_count = 0;
	ambiguous_count = 0; 
	return (mouse_bam, human_bam, graft, host, both, ambiguous, neither, mouse_filtered, human_filtered, total_count,
		host_count, graft_count, both_count, neither_count, ambiguous_count)

def getReads(mouse_filtered, human_filtered, count):
	mouse_read_1 = mouse_filtered[count]
	mouse_read_2 = mouse_filtered[count+1]
	human_read_1 = human_filtered[count]
	human_read_2 = human_filtered[count+1]
	return mouse_read_1, mouse_read_2, human_read_1, human_read_2

def getAlignmentScores(mouse_read_1, mouse_read_2, human_read_1, human_read_2):
	AS_mouse_read_1 = re.findall('(?<=AS:i:)([0-9]*)', mouse_read_1[0])[0]
 	AS_mouse_read_2 = re.findall('(?<=AS:i:)([0-9]*)', mouse_read_2[0])[0]
 	AS_human_read_1 = re.findall('(?<=AS:i:)([0-9]*)', human_read_1[0])[0]
 	AS_human_read_2 = re.findall('(?<=AS:i:)([0-9]*)', human_read_2[0])[0]
 	return AS_mouse_read_1, AS_mouse_read_2, AS_human_read_1, AS_human_read_2

def getNames(mouse_read_1, mouse_read_2, human_read_1, human_read_2):
	mouse_read_1_name = mouse_read_1.split(' ')[0]
	mouse_read_2_name = mouse_read_2.split(' ')[0]
	human_read_1_name = human_read_1.split(' ')[0]
	human_read_2_name = human_read_2.split(' ')[0]

	return mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name

def checkReadNames(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name):
	if mouse_read_1_name == mouse_read_2_name and human_read_1_name == human_read_2_name and mouse_read_1_name == human_read_1_name:
 		sys.exit('read names do not match: M1:{:s} M2:{:s} H1:{:s} H2:{:s}\n'.format(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name))

def classify(AS_mouse_read_1, AS_mouse_read_2, AS_human_read_1, AS_human_read_2, host, host_count, graft, graft_count,
	both, both_count, neither, neither_count, ambiguous, ambiguous_count):
	if (AS_mouse_read_1 > 104 or AS_mouse_read_2 > 104) and AS_human_read_1 < 104 and AS_human_read_2 < 104:
		host.write('{:s}\n'.format(human_read_1_name))
		host_count += 1
	elif AS_mouse_read_1 < 104 and AS_mouse_read_2 < 104 and (AS_human_read_1 > 104 or AS_human_read_2 > 104):
		graft.write('{:s}\n'.format(human_read_1_name))
		graft_count += 1
	elif AS_mouse_read_1 > 104 and AS_mouse_read_2 > 104 and AS_human_read_1 > 104 and AS_human_read_2 > 104:
		both.write('{:s}\n'.format(human_read_1_name))
		both_count += 1
	elif AS_mouse_read_1 < 10 and AS_mouse_read_2 < 10 and AS_human_read_1 < 10 and AS_human_read_2 < 10:
		neither.write('{:s}\n'.format(human_read_1_name))
		neither_count += 1
	else:
		ambiguous.write('{:s}\n'.format(human_read_1_name))
		ambiguous_count += 1
	return graft_count, host_count, both_count, neither_count, ambiguous_count

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
	import sys
	import argparse
	import re

	# argument parser
	# parser = argparse.ArgumentParser(description='Classify reads as host, graft, both, neither, or ambiguous.')
	# parser.add_argument('-m', '--mouse_bam', help='bam file for reads aligned to mouse', type=str, required=True)
	# parser.add_argument('-h', '--human_bam', help='bam file for reads aligned to human', type=str, required=True)
	# parser.add_argument('-o', '--output_dir', help='output directory for output bam and fastq files', type=str, required=True)
	# args = parser.parse_args(sys.argv)
	# mouse_bam = 

	checkArguments()
	(mouse_bam, human_bam, graft, host, both, ambiguous, neither, mouse_filtered, human_filtered, total_count, host_count, 
		graft_count, both_count, neither_count, ambiguous_count) = initializeVariables()
	
	for count in range(0,mouse_filtered,2):

		mouse_read_1, mouse_read_2, human_read_1, human_read_2 = getReads(mouse_filtered, human_filtered, count)
 		AS_mouse_read_1, AS_mouse_read_2, AS_human_read_1, AS_human_read_2 = getAlignmentScores(mouse_read_1, mouse_read_2, human_read_1, human_read_2)
 		mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name = getNames(mouse_read_1, mouse_read_2, human_read_1, human_read_2)
		checkReadNames(mouse_read_1_name, mouse_read_2_name, human_read_1_name, human_read_2_name) 		
        graft_count, host_count, both_count, neither_count, ambiguous_count = classify(AS_mouse_read_1, AS_mouse_read_2, AS_human_read_1, AS_human_read_2, 
        host, host_count, graft, graft_count, both, both_count, neither, neither_count, ambiguous, ambiguous_count)
        total_count += 1

	percentages = calculatePercentages(graft_count, host_count, both_count, ambiguous_count, neither_count, total_count)
	displayOutput(percentages)
		