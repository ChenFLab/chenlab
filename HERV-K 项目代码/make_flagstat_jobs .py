f1 = open('sample.txt', 'r')

namelist = []
for line in f1.readlines():
    line = line.strip().split("/")
    i = line[-1].split('_')[0]
    namelist.append(i)
newnamelist = list(set(namelist))

for name in newnamelist:
    f2 = open(name + 'flag.sh', 'w+')
    f2.write("#!/bin/sh\n#PBS -N "+name+"flag\n#PBS -q core24\n#PBS -o stdout.flag" + name+
    "\n#PBS -e stderr.flag"+ name + "\n#PBS -V\n#PBS -l mem=20gb,walltime=2000:00:00,nodes=1:ppn=8\n#HSCHED -s human+flag+human\n\ncd  /p300s/chenf_group/luhao/20220609-HERV/1-BWA\nsamtools flagstat "+name+".bam > "+name+".flagstat"