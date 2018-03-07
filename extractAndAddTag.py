import sys

array = [line.rstrip('\n') for line in open(sys.argv[1])]
extracted = []
for i in range(len(array)):
    if i % 4 == 0:
        read_name = array[i].split(" ")[0] + '    {:s}'.format(sys.argv[3])
        extracted.append(read_name)
out_file= open(sys.argv[2],'w')
for i in range(len(extracted)):
	out_file.write("{:s}\n".format(extracted[i]))        
