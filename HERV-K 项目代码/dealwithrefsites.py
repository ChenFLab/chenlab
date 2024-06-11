strand_list = open('strand_list.txt','r')
dict1 = {}
for i in strand_list.readlines():
    key = i.strip().split('\t')[0]
    value = i.strip().split('\t')[-1]
    dict1[key] = value
print(dict1)


f1 = open('ref-sites .txt','r')
f2 = open('ref-sites2.txt','a')
for i in f1.readlines():
    start = i.strip().split('\t')[3]
    end = i.strip().split('\t')[4]
    length = int(end)-int(start)
    keyname = i.strip().split('\t')[1]
    strand  = dict1[keyname]
    f2.write(i.strip()+'\t'+str(length)+'\t'+strand+'\n')
