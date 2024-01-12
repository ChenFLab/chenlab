#!/usr/bin/python
import sys
import os
from collections import defaultdict

path = sys.argv[1]
outf = open(sys.argv[2], 'w')
# path = './kegg'
# outf = open('matrix_kegg_raw.csv', 'w')

genome2kegg = {}
genome_list = []
kegg_list = []
files= os.listdir(path)
for file in files: 
     if not os.path.isdir(file): 
          genome = '_'.join(file.split('_')[0:2])
          genome_list.append(genome)
          keggpergenome = defaultdict(int)
          # print(genome)
          f = open(path+"/"+file)
          for line in f.readlines():
               if not line.startswith('#'):                     
                       line1 = line.replace("*"," ").strip(' ').split(' ')[1:]
                       line2 = ' '.join(line1)
                       line = line2.strip(' ').split(' ')
                       # print(line)                      
                       kegg_list.append(line[0])
                       # print(kegg_list)
                       keggpergenome[line[0]] += 1                         
          genome2kegg[genome] = keggpergenome
genome_list = sorted(genome_list)
kegg_list = sorted(list(set(kegg_list)))

header = 'genes,' + ','.join(genome_list) + '\n'
outf.write(header)
for kegg in kegg_list:
     line = []
     for genome in genome_list:
          # print(genome)
          keggpergenome = genome2kegg[genome]
          if kegg in keggpergenome.keys():
               line.append(str(keggpergenome[kegg]))
          else:
               line.append('0')
     line = kegg + ',' + ','.join(line) + '\n'
     outf.write(line)
