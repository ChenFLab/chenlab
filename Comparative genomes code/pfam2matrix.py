#!/usr/bin/python
import sys
import os
from collections import defaultdict

path = sys.argv[1]
outf = open(sys.argv[2], 'w')
# path = './pfam_test'
# outf = open('matrix_pfam_raw.csv', 'w')

genome2pfam = {}
genome_list = []
pfam_list = []
files= os.listdir(path)
for file in files: 
     if not os.path.isdir(file): 
          genome = '_'.join(file.split('_')[0:2])
          genome_list.append(genome)
          pfampergenome = defaultdict(int)
          # print(genome)
          f = open(path+"/"+file)
          for line in f.readlines():
               if not line.startswith('#') and line != '\n': 
                    line = line.strip().split(' ')
                    line = [value for value in line if value != '']
                    if len(line) != 15:
                         print(line)
                    else:
                         pfampergenome[line[5]] += 1
                         pfam_list.append(line[5])
          genome2pfam[genome] = pfampergenome
# print(genome2pfam)

genome_list = sorted(genome_list)
pfam_list = sorted(list(set(pfam_list)))

header = 'genes,' + ','.join(genome_list) + '\n'
outf.write(header)

for pfam in pfam_list:
     line = []
     for genome in genome_list:
          # print(genome)
          pfampergenome = genome2pfam[genome]
          if pfam in pfampergenome.keys():
               line.append(str(pfampergenome[pfam]))
          else:
               line.append('0')
     line = pfam + ',' + ','.join(line) + '\n'
     outf.write(line)
