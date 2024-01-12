#!/usr/bin/python
import sys
from collections import defaultdict

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
outf = open(sys.argv[3], 'w')

# inf1 = open('all_info_prokka.txt', 'r')
# inf2 = open('Orthogroups.txt', 'r')
# outf = open('orthogroup_info_prokka.txt', 'w')

info_dict = {}
for line in inf1.readlines():
	line = line.strip('\n').split('\t')
	info_dict[line[0]] = [line[3], line[6], line[5], line[4]]

for line in inf2.readlines():
	gene_protein_list = []
	COG_list = []
	EC_list = []
	line = line.strip('\n').split(' ')
	orthogroup = line[0].strip(':')
	acc_list = line[1:]
	for acc in acc_list:
		if acc in info_dict:
			if info_dict[acc][0] != '':
				gene_protein_list.append(info_dict[acc][0] + '\t' + info_dict[acc][1])
			if info_dict[acc][2] != '':
				COG_list.append(info_dict[acc][2])
			if info_dict[acc][3] != '':
				EC_list.append(info_dict[acc][3])
	gene_protein_list = sorted(list(set(gene_protein_list)))
	gene_list = []
	protein_list = []
	for item in gene_protein_list:
		item = item.split('\t')
		gene_list.append(item[0])
		protein_list.append(item[1])
	gene = ' | '.join(gene_list)
	protein = ' | '.join(protein_list)
	COG = ' | '.join(sorted(list(set(COG_list))))
	EC = ' | '.join(sorted(list(set(EC_list))))
	outf.write(orthogroup + '\t' + gene + '\t' + protein + '\t' + COG + '\t' + EC + '\n')

inf1.close()
inf2.close()
outf.close()
