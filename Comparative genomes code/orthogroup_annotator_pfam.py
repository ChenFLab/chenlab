#!/usr/bin/python
import sys
from collections import defaultdict

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
inf3 = open(sys.argv[3], 'r')
outf = open(sys.argv[4], 'w')


# inf1 = open('Pfam-A.clans.tsv', 'r')
# inf2 = open('all_info_pfam.txt', 'r')
# inf3 = open('Orthogroups.txt', 'r')
# outf = open('orthogroup_info_pfam.txt', 'w')


pfam_dict = {}
for line in inf1.readlines():
	line = line.strip().split('\t')
	info = [line[3], line[4], line[1], line[2]]
	pfam_dict[line[0]] = info
# print(pfam_dict)


conversion_dict = defaultdict(list)
for line in inf2.readlines():
	line = line.strip().split(' ')
	line = [value for value in line if value != '']
	conversion_dict[line[0]].append(line[5].split('.')[0])
# print(conversion_dict)


for line in inf3.readlines():
	line = line.strip('\n').split(' ')
	orthogroup = line[0].strip(':')
	acc_list = line[1:]
	pfam_list = []
	for acc in acc_list:
		if acc in conversion_dict:
			pfam = sorted(conversion_dict[acc])
			for pfam_acc in pfam:
				pfam_list.append(pfam_acc)
	pfam_list = sorted(list(set(pfam_list)))
	# print(pfam_list)
	l1 = []
	l2 = []
	l3 = []
	l4 = []
	for pfam_acc in pfam_list:
		l1.append(pfam_dict[pfam_acc][0])
		l2.append(pfam_dict[pfam_acc][1])
		l3.append(pfam_dict[pfam_acc][2])
		l4.append(pfam_dict[pfam_acc][3])
	pfam_s = ' | '.join(pfam_list)
	s1 = ' | '.join(l1)
	s2 = ' | '.join(l2)
	if set(l3) != {''}:
		s3 = ' | '.join(l3)
		s4 = ' | '.join(l4)
	else:
		s3 = ''
		s4 = ''
	outf.write(orthogroup + '\t' + pfam_s + '\t' + 
		s1 + '\t' + s2 + '\t' + s3 + '\t' + s4 + '\n')


inf1.close()
inf2.close()
inf3.close()
outf.close()
