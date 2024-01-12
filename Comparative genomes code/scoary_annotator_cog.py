#!/usr/bin/python
import sys

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
inf3 = open(sys.argv[3], 'r')
outf = open(sys.argv[4], 'w')

# inf1 = open('fun-20.rev.tab', 'r')
# inf2 = open('cog-20.def.rev.tab', 'r')
# inf3 = open('scoary_cog_2763.results.csv', 'r')
# outf = open('scoary_cog_2763.results.ann.csv', 'w')


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
	# l2 = list(set(l2))
	s1 = ' | '.join(l1)
	s2 = ' | '.join(l2)
	info = [line[3], line[2], line[4], line[1], s1, s2]
	for i in range(len(info)):
		info[i] = '"' + info[i] + '"'
	cog_dict[line[0]] = ','.join(info)
print(cog_dict)

header = inf3.readline()
header = header.strip('\n')
outf.write(header + ',Gene_symbol,Gene_description,Gene_function,COG_class,COG_class_description,COG_class_summary\n')
for line in inf3.readlines():
	line = line.strip().split(',')
	acc = line[0].replace('"', '').split('.')[0]
	# print(acc)
	# print(cog_dict[acc])
	line = ','.join(line) + ',' + cog_dict[acc] + '\n'
	outf.write(line)

inf1.close()
inf2.close()
inf3.close()
outf.close()
