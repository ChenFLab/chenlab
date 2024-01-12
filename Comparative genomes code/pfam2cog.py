#!/usr/bin/python
import sys
import os
from collections import defaultdict

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
outf = open(sys.argv[3], 'w')
# path = './pfam_test'
# outf = open('matrix_pfam_raw.csv', 'w')



gene_pfam_dict = defaultdict(list)

for line in inf1.readlines(): 
    if not line.startswith('#') and line != '\n': 
        line = line.strip().split(' ')
        line = [value for value in line if value != '']
        if len(line) != 15:
            print(line)
        else:
            gene_name = line[0]
            pfam_id = line[5]
            gene_pfam_dict[gene_name].append(pfam_id)
        #print(gene_dict)


gene_cog_dict = {}
protein_list = []

for line in inf2.readlines():
    line = line.strip().split('\t')
    if line[0] not in protein_list:
        if abs((float(line[7]) - float(line[6])) / float(line[8])) >= 0.7:
            protein_list.append(line[0])
            protein_name =line[0]
            cog_id = line[-1].split(', ')[0]
            gene_cog_dict[protein_name] = cog_id


gene_dict = defaultdict(list)
for gene_name in gene_pfam_dict.keys():
    gene_dict[gene_name] = gene_pfam_dict[gene_name]
    if gene_name in gene_cog_dict:
        gene_dict[gene_name].append(gene_cog_dict[gene_name])
    else:
        gene_dict[gene_name].append("NA")


for gene_name, values in gene_dict.items():
    outf.write(f"{gene_name}\t{'|'.join(values)}\n")

inf1.close()
inf2.close()
outf.close()