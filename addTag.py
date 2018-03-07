import sys

array = [line.rstrip('\n') for line in open(sys.argv[1])]
for i in len(array):
	array[i] = array[i] + '    {:s}'.format(sys.argv[2])
out_file= open('tagged_{:s}.fastq'.format(sys.argv[2]),'w')
for i in len(array):
	out_file.write("{:s}\n".format(array[i]))	
