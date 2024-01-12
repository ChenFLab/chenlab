#!/usr/bin/python
import sys
from collections import defaultdict

inf = open(sys.argv[1], 'r')
outf1 = open(sys.argv[2], 'w')
outf2 = open(sys.argv[3], 'w')
threshld = sys.argv[4]


# inf = open('distances.txt', 'r')
# outf1 = open('distances_clean.txt', 'w')
# outf2 = open('distances_filtered.txt', 'w')

for line in inf.readlines():
	line = line.strip().split('\t')
	
	strain1 = line[0].split('/')[-1].replace('.fna', '')
	strain2 = line[1].split('/')[-1].replace('.fna', '')
	if strain1 != strain2: 
		combo = '\t'.join(sorted([strain1,strain2]))
		outf1.write(combo + '\t' + '\t'.join(line[2:]) + '\n')
		if float(line[2]) <= float(threshld):		
			outf2.write(combo +'\t' + line[2] + '\n')

inf.close()
outf1.close()
outf2.close()