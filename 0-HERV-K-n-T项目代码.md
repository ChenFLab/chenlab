1 bwa生成bam文件

```shell
bwa mem -t 4 -R '@RG\tID:101-5\tPL:illumina\tLB:101-5\tPU:101-5\tSM:101-5' /p300s/chenf_group/luhao/20220329-HERV/idex/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /asnas/chenf_group/licd/HERVdata/20220315/2.cleandata/101-5_FDHG220033325-1a/101-5_FDHG220033325-1a_1.clean.fq.gz /asnas/chenf_group/licd/HERVdata/20220315/2.cleandata/101-5_FDHG220033325-1a/101-5_FDHG220033325-1a_2.clean.fq.gz | samtools view -Sb - > /p300s/chenf_group/luhao/20220609-HERV/1-BWA//101-5.bam
```

2 排序 并计算测序深度

```shell 
cd  /p300s/chenf_group/luhao/20220609-HERV/1-BWA
samtools sort -@ 4 -m 4G -o 101-5_sorted.bam 101-5.bam

# 计算测序深度
samtools depth 101-5_sorted.bam > 101-5_sorted.depth
```

3  使用bed文件根据对应坐标提取mapped reads 

```shell
cd  /p300s/chenf_group/luhao/20220609-HERV/1-BWA
samtools view -b -L /p300s/chenf_group/luhao/20220329-HERV/polymorphicHERV/polymorphicHERV-master/bed_files/HERVK_hg38_sort_01apr2018.bed 101-5_sorted.bam > 101-5_sorted_extracted.bam
```

4 将bam文件转换为fasta文件 

```shell
cd /p300s/chenf_group/luhao/20220609-HERV/1-BWA
samtools view 101-5_sorted_extracted.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' > 101-5_sorted_extracted.fa
```

5 使用dsk软件计数，  是一种快速、高效地计算 DNA 中 k-mer 的工具。k-mer 是长度为 k 的子串，计算所有此类子串的出现次数是许多 DNA 序列分析的核心步骤。 

```shell
cd /p300s/chenf_group/luhao/20220329-HERV/polymorphicHERV/polymorphicHERV-master
/p300s/chenf_group/luhao/software/dsk-v2.2.0-bin-Linux/bin/dsk -file 101-5_sorted_extracted.fa -kmer-size 50 -abundance-min 0
/p300s/chenf_group/luhao/software/dsk-v2.2.0-bin-Linux/bin/dsk2ascii -file 101-5__sortedextracted.h5 -out 101-5_sorted_extracted.txt
perl tofasta.pl /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5_sorted_extracted.txt /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5_sorted_extracted.50.fa 
```

6  精确匹配数据到参考，然后生成n/T矩阵

```shell
cd /p300s/chenf_group/luhao/20220329-HERV/polymorphicHERV/polymorphicHERV-master
python3 findmatch.py unique.withrc.50.fa /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5_sorted_extracted.50.fa /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5.dat
python3 labelcount.py /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5.dat sortedSites /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5.label
python3 concate_matrics.py T.50 /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5.label /p300s/chenf_group/luhao/20220609-HERV/1-BWA/101-5.50.dat
```





