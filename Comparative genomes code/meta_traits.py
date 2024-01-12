import sys
import os
from collections import defaultdict

inf1 = open(sys.argv[1], 'r')
outf1 = open(sys.argv[2], 'w')

genome_list = []
for line in inf1.readlines():      
     line = ','.join(line.split('\t'))
     line = line.strip('\n')
     genome_list.append(line)
print(genome_list)
outf1.write(',TraitY' + '\n' +'\n'.join(genome_list))