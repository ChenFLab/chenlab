#!/usr/bin/python
import sys

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
outf1 = open(sys.argv[3], 'w')
outf2 = open(sys.argv[4], 'w')

# inf1 = open('traits_growth_test.csv', 'r')
# inf2 = open('matrix_pfam_raw.csv', 'r')
# outf1 = open('matrix_pfam_growth_numeric.csv', 'w')
# outf2 = open('matrix_pfam_growth_binary.csv', 'w')

genome_list = []
inf1.readline()
for line in inf1.readlines():
	line = line.strip().split(',')
	genome_list.append(line[0])
# print(genome_list)
# print(len(genome_list))

outf1.write('genes,' + ','.join(genome_list) + '\n')
outf2.write('genes,' + ','.join(genome_list) + '\n')

count_dict = {}
header = inf2.readline()
header = header.strip().split(',')
header = header[1:]
# print(len(header))

for line in inf2.readlines():
	line = line.strip().split(',')
	protein = line[0]
	count_list = line[1:]
	numeric_line = []
	binary_line = []
	for genome in genome_list:
		ind = header.index(genome)
		numeric_val = count_list[ind]
		if int(numeric_val) == 0:
			binary_val = 0
		elif int(numeric_val) > 0:
			binary_val = 1
		else:
			print('ERROR')
		numeric_line.append(numeric_val)
		binary_line.append(str(binary_val))
	# print(numeric_line)
	# print(binary_line)
	numeric_line = protein + ',' + ','.join(numeric_line) + '\n'
	binary_line = protein + ',' + ','.join(binary_line) + '\n'
	# print(numeric_line)
	# print(binary_line)
	outf1.write(numeric_line)
	outf2.write(binary_line)

inf1.close()
inf2.close()
outf1.close()
outf2.close()
