f1 = open('sample.txt', 'r')

namelist = []
for line in f1.readlines():
    line = line.strip().split("/")
    i = line[-1].split('_')[0]
    namelist.append(i)
print(namelist)
newnamelist = list(set(namelist))
print(len(newnamelist))

for name in newnamelist:
    f2 = open(name + 'whole_procedure.sh', 'w+')
    f2.write("#!/bin/sh\n#PBS -N "+name+"whole procedure\n#PBS -q core24\n#PBS -o stdout." + name+
    "\n#PBS -e stderr."+ name + "\n#PBS -V\n#PBS -l mem=20gb,walltime=2000:00:00,nodes=1:ppn=8\n#HSCHED -s human+nT+human\n\ncd  /p300s/chenf_group/luhao/20220609-HERV/1-BWA\nsamtools sort -@ 4 -m 4G -o " + name +"_sorted.bam " +
    name + ".bam\nsamtools view -b -L /p300s/chenf_group/luhao/20220329-HERV/polymorphicHERV/polymorphicHERV-master/bed_files/HERVK_hg38_sort_01apr2018.bed " + name +"_sorted.bam > "+name+
    "_sorted_extracted.bam\nsamtools view "+name+r'''_sorted_extracted.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' > '''+name+"_sorted_extracted.fa\n/p300s/chenf_group/luhao/software/dsk-v2.2.0-bin-Linux/bin/dsk -file "+name+"_sorted_extracted.fa -kmer-size 50 -abundance-min 0\n/p300s/chenf_group/luhao/software/dsk-v2.2.0-bin-Linux/bin/dsk2ascii -file "+name+"_sorted_extracted.h5 -out "+name+"_sorted_extracted.txt\ncd /p300s/chenf_group/luhao/20220329-HERV/polymorphicHERV/polymorphicHERV-master\nperl tofasta.pl /p300s/chenf_group/luhao/20220609-HERV/1-BWA/"+name+"_sorted_extracted.txt /p300s/chenf_group/luhao/20220609-HERV/1-BWA/"+name+"_sorted_extracted.50.fa\npython3 findmatch.py unique.withrc.50.fa /p300s/chenf_group/luhao/20220609-HERV/1-BWA/"+name+"_sorted_extracted.50.fa /p300s/chenf_group/luhao/20220609-HERV/1-BWA/"+name+'.dat\npython3 labelcount.py /p300s/chenf_group/luhao/20220609-HERV/1-BWA/'+name+".dat sortedSites /p300s/chenf_group/luhao/20220609-HERV/1-BWA/"+name+'.label\npython3 concate_matrics.py T.50 /p300s/chenf_group/luhao/20220609-HERV/1-BWA/'+name+".label /p300s/chenf_group/luhao/20220609-HERV/1-BWA/"+name+'.50.dat\n')
