"""
extractAS.py
Created by Heather D'Souza

Write AS score output in an R-friendly format.
Created by Heather D'Souza
"""
import sys
import re

if len(sys.argv) < 3:
	print("Usage: python extractAS.py diff_class.txt class_AS.txt")
	sys.exit()

in_file = open(sys.argv[2], 'r')
alignment_file = open(sys.argv[3], 'w')
with open(sys.argv[1], 'r') as in_file:
	for line in in_file:
		as_tags = re.findall('AS:i:\S+', line, flags=0)
		as_numbers = []
		for i in as_tags:
			as_numbers.append(i.split(':')[2])
		alignment_file.write("{:3s} {:3s}".format(as_numbers[0], as_numbers[1]))
alignment_file.close()
in_file.close()
	