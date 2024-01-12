#!/usr/bin/python
import sys
from collections import defaultdict

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
inf3 = open(sys.argv[3], 'r')
inf4 = open(sys.argv[4], 'r')
outf = open(sys.argv[5], 'w')


# inf1 = open('fun-20.rev.tab', 'r')
# inf2 = open('cog-20.def.rev.tab', 'r')
# inf3 = open('all_info_cog.txt', 'r')
# inf4 = open('Orthogroups.txt', 'r')
# outf = open('orthogroup_info_cog.txt', 'w')

fun_dict = {}
for line in inf1.readlines():
	line = line.strip('\n').split('\t')
	# print(line)
	line[2] = line[2].replace('"', '')
	info = [line[2], line[3]]
	fun_dict[line[0]] = info
# print(fun_dict)

cog_dict = {}
for line in inf2.readlines():
	line = line.strip('\n').replace('"', '').split('\t')
	items = line[1]
	l1 = []
	l2 = []
	for item in items:
		l1.append(fun_dict[item][0])
		l2.append(fun_dict[item][1])
	s1 = '; '.join(l1)
	s2 = '; '.join(l2)
	info = [line[3], line[2], line[4], line[1], s1, s2]
	cog_dict[line[0]] = info
# print(cog_dict)

conversion_dict = {}
flag = ''
for line in inf3.readlines():
	line = line.strip('\n').split('\t')
	if line[0] != flag:
		if abs((float(line[7]) - float(line[6])) / float(line[8])) >= 0.7:
			conversion_dict[line[0]] = line[11].split(',')[0]
			flag = line[0]
# print(conversion_dict)

for line in inf4.readlines():
	line = line.strip('\n').split(' ')
	orthogroup = line[0].strip(':')
	acc_list = line[1:]
	cog_list = []
	for acc in acc_list:
		if acc in conversion_dict:
			cog = conversion_dict[acc]
			cog_list.append(cog)
	cog_list = sorted(list(set(cog_list)))
	# print(cog_list)
	l1 = []
	l2 = []
	l3 = []
	l4 = []
	l5 = []
	l6 = []
	for cog_acc in cog_list:
		l1.append(cog_dict[cog_acc][0])
		l2.append(cog_dict[cog_acc][1])
		l3.append(cog_dict[cog_acc][2])
		l4.append(cog_dict[cog_acc][3])
		l5.append(cog_dict[cog_acc][4])
		l6.append(cog_dict[cog_acc][5])
	cog_s = ' | '.join(cog_list)
	s1 = ' | '.join(l1)
	s2 = ' | '.join(l2)
	s3 = ' | '.join(l3)
	s4 = ' | '.join(l4)
	s5 = ' | '.join(l5)
	s6 = ' | '.join(l6)
	outf.write(orthogroup + '\t' + cog_s + '\t' + 
		s1 + '\t' + s2 + '\t' + s3 + '\t' + s4 + '\t' + s5 + '\t' + s6 + '\n')



inf1.close()
inf2.close()
inf3.close()
inf4.close()
outf.close()
