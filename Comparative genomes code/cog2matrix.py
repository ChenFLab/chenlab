#!/usr/bin/python
import sys
import os
from collections import defaultdict

path = sys.argv[1]
outf = open(sys.argv[2], 'w')
# path = './cog_test'
# outf = open('matrix_cog_raw.csv', 'w')

genome2cog = {}
genome_list = []
cog_list = []
files= os.listdir(path)
for file in files: 
     if not os.path.isdir(file): 
          genome = '_'.join(file.split('_')[0:2])
          genome_list.append(genome)
          cogpergenome = defaultdict(int)
          protein_list = []
          # print(genome)
          f = open(path+"/"+file)
          for line in f.readlines():
               line = line.strip().split('\t')
               if line[0] not in protein_list:
                    if abs((float(line[7]) - float(line[6])) / float(line[8])) >= 0.7:
                         protein_list.append(line[0])
                         cog = line[-1].split(', ')[0]
                         cogpergenome[cog] += 1
                         cog_list.append(cog)
          genome2cog[genome] = cogpergenome
# print(genome2cog)

genome_list = sorted(genome_list)
cog_list = sorted(list(set(cog_list)))

header = 'genes,' + ','.join(genome_list) + '\n'
outf.write(header)

for cog in cog_list:
     line = []
     for genome in genome_list:
          # print(genome)
          cogpergenome = genome2cog[genome]
          if cog in cogpergenome.keys():
               line.append(str(cogpergenome[cog]))
          else:
               line.append('0')
     line = cog + ',' + ','.join(line) + '\n'
     outf.write(line)
