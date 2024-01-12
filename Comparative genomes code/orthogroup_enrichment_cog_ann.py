#!/usr/bin/python
import sys
import os
from collections import defaultdict

inf1 = open(sys.argv[1], 'r')
inf2 = open(sys.argv[2], 'r')
outf = open(sys.argv[3], 'w')

# inf1 = '/workdir/zhangml/20221024_NTM/restart/13.orthpgroup_enrichment/1.1.1.cog/OG0007945.cog'
# inf2 = '/workdir/zhangml/20221024_NTM/restart/13.orthpgroup_enrichment/1.1.1.kegg/OG0007945.kegg'
# outf = open('/workdir/zhangml/20221024_NTM/restart/13.orthpgroup_enrichment/ann/1.1.1.ann/OG0007945.ann.csv', 'w')


cog_list = []
cogpergenome = defaultdict(int)
protein_list = []
for line in inf1.readlines():
       line = line.strip().split('\t')
       #print(line)
       if line[0] not in protein_list:
            if abs((float(line[7]) - float(line[6])) / float(line[8])) >= 0.7:
                 protein=line[0]
                 protein_list.append(protein)
                 cog = line[-1]
                 cog_list.append(cog)
                 cogpergenome[protein] = cog
keggpergenome = defaultdict(int)
kegg_list = []
protein_list = []
for line in inf2.readlines():
       if not line.startswith('#'):
           line = line.replace("*","\t")
           line=line.strip('\t').strip('\n').split('\t')
           #print(line[0])
           if line[0] not in protein_list:
               protein=line[0]
               #print(protein)
               protein_list.append(protein)
               #print(protein_list)
               kegg=line[-1]
               kegg_list.append(kegg)
               keggpergenome[protein] = kegg


header = 'genes,' + 'cog,'+  'kegg,' + '\n'
outf.write(header)
line =[]
for protein in protein_list:
       cog = cogpergenome[protein]
       kegg = keggpergenome[protein]
       line = protein + ',' + str(cog)+ ','+str(kegg)+ ',' + '\n'
       outf.write(line)
inf1.close()
inf2.close()
outf.close()