#!/usr/bin/python
import sys

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
outf = open(sys.argv[3], 'w')

# inf1 = open('Pfam-A.clans.tsv', 'r')
# inf2 = open('scoary_pfam_2763.results.csv', 'r')
# outf = open('scoary_pfam_2763.results.ann.csv', 'w')

pfam_dict = {}
for line in inf1.readlines():
	line = line.strip().split('\t')
	info = [line[3], line[4], line[1], line[2]]
	for i in range(len(info)):
		info[i] = '"' + info[i] + '"'
	pfam_dict[line[0]] = ','.join(info)
# print(pfam_dict)

header = inf2.readline()
header = header.strip('\n')
outf.write(header + ',Gene_symbol,Gene_description,Pfam_clan,Pfam_clan_description\n')
for line in inf2.readlines():
	line = line.strip().split(',')
	acc = line[0].replace('"', '').split('.')[0]
	# print(acc)
	# print(pfam_dict[acc])
	line = ','.join(line) + ',' + pfam_dict[acc] + '\n'
	outf.write(line)

inf1.close()
inf2.close()
outf.close()
