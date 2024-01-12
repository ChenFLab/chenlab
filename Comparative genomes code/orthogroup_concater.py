#!/usr/bin/python
import sys
import os
from collections import defaultdict

inf = open(sys.argv[1], 'r')
path = sys.argv[2]
outf = open(sys.argv[3], 'w')

# inf = open('prokka_table_366.txt', 'r')
# path = './alignments'
# outf = open('orthologues_singlecopy.aln', 'w')


prokka_dict = {}
for line in inf.readlines():
     line = line.strip('\n').split('\t')
     # print(line)
     prokka_dict[line[1]] = line[0]
inf.close()

concatenate_aln_dict = defaultdict(str)
files= os.listdir(path)
for file in files: 
     if not os.path.isdir(file): 
          f = open(path+"/"+file)
          for line in f.readlines():
               line = line.strip('\n')
               if line.startswith('>'):
                    name = line.strip('>').split('_')[0]
                    # print(name)
                    name = prokka_dict[name]
               else:                
                    concatenate_aln_dict[name] += line
# print(concatenate_aln_dict)

for key in sorted(concatenate_aln_dict):
     outf.write('>' + key + '\n' + concatenate_aln_dict[key] + '\n')
outf.close() 