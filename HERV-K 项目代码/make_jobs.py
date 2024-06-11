all_file = open("114-sample-md5.txt","r")


output_path = "/p300s/chenf_group/luhao/20220609-HERV/1-BWA/"

sample = []

for line in all_file:
    inf = line.strip().split("  ")
    path = inf[1]
    if "clean" in inf[1]:
        file_name = path.split("/")[-1]
        ID = file_name.split("_")[0]
        if ID not in sample:
            sample.append(ID)
            tmp_out = open(output_path+ID+"_bwa_job.sh","w")
            tmp_out.write("#!/bin/sh\n#PBS -N bwa_"+ID+"\n#PBS -q core24\n"+"#PBS -o stdout."+ID.replace("-","_")+
                          "\n#PBS -e stderr."+ID.replace("-","_")+
                          "\n#PBS -V\n#PBS -l mem=20gb,walltime=2000:00:00,nodes=1:ppn=7\n#HSCHED -s human+samtools+human\ncd "+
                          output_path+"\n"+
                          "bwa mem -t 4 -R '@RG\\tID:"+ID+"\\tPL:illumina\\tLB:"+ID+"\\tPU:"+ID+"\\tSM:"+ID+"' "+
                          "/p300s/chenf_group/luhao/20220329-HERV/idex/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna "+
                          inf[1]+" "+inf[1].replace("1.clean.fq.gz","2.clean.fq.gz")+" | samtools view -Sb - > "+output_path+"/"+ID+".bam")
            tmp_out.close()
            
all_file.close()
