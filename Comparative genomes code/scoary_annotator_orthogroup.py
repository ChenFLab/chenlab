#!/usr/bin/python
import sys

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
inf3 = open(sys.argv[3], 'r')
inf4 = open(sys.argv[4], 'r')
outf = open(sys.argv[5], 'w')

# inf1 = open('orthogroup_info_prokka.txt', 'r')
# inf2 = open('orthogroup_info_cog.txt', 'r')
# inf3 = open('orthogroup_info_pfam.txt', 'r')
# inf4 = open('scoary_orthogroup_366_rev.results.csv', 'r')
# outf = open('scoary_orthogroup_366_rev.results.ann.csv', 'w')

prokka_dict = {}
for line in inf1.readlines():
	line = line.strip('\n').split('\t')
	info = line[1:3]
	for i in range(len(info)):
		info[i] = '"' + info[i] + '"'
	prokka_dict[line[0]] = ','.join(info)
# print(prokka_dict)

cog_dict = {}
for line in inf2.readlines():
	line = line.strip('\n').split('\t')
	info = line[1:]
	for i in range(len(info)):
		info[i] = '"' + info[i] + '"'
	cog_dict[line[0]] = ','.join(info)
# print(cog_dict)

pfam_dict = {}
for line in inf3.readlines():
	line = line.strip('\n').split('\t')
	info = line[1:]
	for i in range(len(info)):
		info[i] = '"' + info[i] + '"'
	pfam_dict[line[0]] = ','.join(info)
# print(pfam_dict)

header = inf4.readline()
header = header.strip('\n')
outf.write(header + ',Prokka_gene_symbol,Prokka_protein_product,\
	COG_id,COG_gene_symbol,COG_gene_description,COG_gene_function,\
	COG_class,COG_class_description,COG_class_summary,\
	Pfam_id,Pfam_domain_symbol,Pfam_domain_description,Pfam_clan,Pfam_clan_description\n')
for line in inf4.readlines():
	line = line.strip().split(',')
	# print(line)
	acc = line[0].replace('"', '')
	# print(prokka_dict[acc])
	line = ','.join(line) + ',' + prokka_dict[acc] + ',' + cog_dict[acc] + ',' + pfam_dict[acc] + '\n'
	outf.write(line)

inf1.close()
inf2.close()
inf3.close()
inf4.close()
outf.close()
