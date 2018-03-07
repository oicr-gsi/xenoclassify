"""
extractReadName.py
Created by Heather D'Souza
"""
import sys

def get_name(infile):
    """ extract read names from fastq """
    o_file = sys.argv[2]
    read_names = open(o_file, 'w')
    count = 0
    with open(infile) as fastq:
        for line in fastq:
            count += 1
            if count % 4 == 1:
                line = line.split(" ")[0]
                read_names.write(line + '\n')

if __name__ == "__main__":
    get_name(sys.argv[1])
