"""
diff.py
Created by Heather D'Souza

GNU diff tries to minimizes the number of differences between files, by
shifting down lines of a file if it believes there is a "gap" in the other 
file.

This diff script will compare files line by line and output the lines in each file side-by-side.
If there is a difference (disregarding white spaces and new lines), a "|" will be inserted between the two lines.

This function assumes that files have the same number of newlines.
"""
import sys

if len(sys.argv) < 3:
	print("Usage: python diff.py file1 file2")
	sys.exit()

file1 = open(sys.argv[1],'r')
file2 = open(sys.argv[2],'r')

for line1, line2 in zip([line.strip() for line in file1], [line.strip() for line in file2]):	
	if line1 != line2:
		sys.stdout.write("{:70s} {:3s} {:70s}\n".format(line1, '|', line2)) 
	else:
		sys.stdout.write("{:70s} {:3s} {:70s}\n".format(line1, ' ', line2))  
		