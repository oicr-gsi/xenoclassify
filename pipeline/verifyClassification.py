"""
verifyClassification.py
Created by Heather D'Souza

This script determines how similarly reads are classified by our classifier and Xenome.
"""

def checkArguments ():
	if len(sys.argv) != 6:
		print "Usage: python verifyClassification.py xenome.txt graft.txt host.txt both.txt neither.txt"

def indexFile (in_file):
	indexed_file = [line.rstrip('\n') for line in open(in_file)]
	file_length = len(indexed_file)
	return indexed_file, file_length

def verifyClass (class_counter, tag, class_name, correct, incorrect, class_correct):
	class_counter += 1
	if tag==class_name:
		correct += 1
		class_correct += 1
	else:
		incorrect += 1
	return class_counter, correct, incorrect, class_correct

def countReads (totals, tag):
	if tag == 'graft':
		totals[0] += 1
	elif tag == 'host':
		totals[1] += 1
	elif tag == 'both':
		totals[2] += 1
	elif tag == 'ambiguous':
		totals[3] += 1
	elif tag == 'neither':
		totals[4] += 1
	return totals

# checks if counter is greater than array length
def checkCounter (counter, length):
	if counter > length - 1:
		counter = 0
	return counter

def outputResults (correct, incorrect, totals, g_correct, h_correct, b_correct, a_correct, n_correct):
	total=correct+incorrect
	percentage_correct=round(float(correct)/total*100,1)
	percentage_incorrect=round(float(incorrect)/total*100,1)
	g_correct=round(float(g_correct)/totals[0]*100,1)
	h_correct=round(float(h_correct)/totals[1]*100,1)
	b_correct=round(float(b_correct)/totals[2]*100,1)
	n_correct=round(float(n_correct)/totals[4]*100,1)
	print "Percentage of Reads Classified Correctly\n"
	print "Total: {:.1f}%".format(percentage_correct)
	print "Graft: {:.1f}%".format(g_correct)
	print "Host: {:.1f}%".format(h_correct)
	print "Both: {:.1f}%".format(b_correct)
	print "Neither: {:.1f}%".format(n_correct)

if __name__ == '__main__':
	
	import re
	import sys 	

	checkArguments()	

	xenome = open(sys.argv[1],'r')

	# create array for input files
	graft, g_length = indexFile(sys.argv[2])
	host, h_length = indexFile(sys.argv[3])
	both, b_length = indexFile(sys.argv[4])
	ambiguous, a_length = ([], 0)
	neither, n_length = indexFile(sys.argv[5])

	# counters for arrays
	g_count = 0
	h_count = 0
	b_count = 0
	a_count = 0
	n_count = 0
	correct = 0
	incorrect = 0
	g_correct = 0
	h_correct = 0
	b_correct = 0
	a_correct = 0
	n_correct = 0

	# counter for length of files
	# number of reads in Xenome's classes
	totals = [0,0,0,0,0]

	# traverse through xenome read names
	for line in xenome:

		g_count = checkCounter(g_count,g_length)
		h_count = checkCounter(h_count,h_length)
		b_count = checkCounter(b_count,b_length)
		a_count = 0
		n_count = checkCounter(n_count,n_length)

		# separate read name and tag
		line = line.strip()
		contents = line.split()
		read_name = contents[0]
		tag = contents[1]

		# count total reads in Xenome class
		totals = countReads(totals, tag)

		# verify classificiation
		if read_name == graft[g_count]:
			g_count, correct, incorrect, g_correct = verifyClass(g_count, tag, "graft", correct, incorrect, g_correct)
		elif read_name == host[h_count]:
			h_count, correct, incorrect, h_correct = verifyClass(h_count, tag, "host", correct, incorrect, h_correct)
		elif read_name == both[b_count]:
			b_count, correct, incorrect, b_correct = verifyClass(b_count, tag, "both", correct, incorrect, b_correct)
		elif read_name == neither[n_count]:
			n_count, correct, incorrect, n_correct = verifyClass(n_count, tag, "neither", correct, incorrect, n_correct)
		else:
			incorrect += 1

	outputResults(correct, incorrect, totals, g_correct, h_correct, b_correct, a_correct, n_correct)
