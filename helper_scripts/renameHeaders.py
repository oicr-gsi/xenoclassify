"""
renameHeaders.py
Created by Heather D'Souza

Before a hybrid index can be created from hg19_random.fa and 
mm10.fa, the headers must be changed so that they are reference-specific.
This is important during the alignment process, as the name of the chromosome to which
each read aligns is included in the SAM format for that read. The only way to
distinguish if the read aligned to i.e. chr1 on mouse or chr1 on human is by the
reference name in the chromosome tag.
"""
import re
import sys 

out_file = open('hg19_random_tagged.fa', 'w')

with open(sys.argv[1], 'r') as in_file:
	for line in in_file:
		if line.strip():
			chr_name = re.findall('>chr[0-9]*', line, flags=0)
			if len(chr_name) > 0:
				chr_name = chr_name + '_hg19'
				out_file.write("{:s}\n".format(chr_name))
				continue
		out_file.write("{:s}\n".format(line))
