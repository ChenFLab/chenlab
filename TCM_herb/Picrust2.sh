#Picrust2 功能预测
wd=/workdir/wangxt/EasyAmpliconmaster/picrust2
# 运行流程
    picrust2_pipeline.py -s otus.fa -i otutab.txt \
      -o picrust2 -p 8
    # 添加EC/KO/Pathway注释
    add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
    add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
      -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 
  
# METACYC pathway的注释 
 add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
      -o pathways_out/path_abun_unstrat_descrip.tsv.gz  

# KO映射KEGG pathway 
pathway_pipeline.py -i picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -o picrust2/KEGG_pathways --no_regroup --map ~/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv


#EC 映射METACYC pathway
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -o pathways_out \
                    --intermediate minpath_working \
                    -p 1
# https://github.com/picrust/picrust2/wiki/Infer-pathway-abundances

db=/workdir/wangxt/db
script=/workdir/wangxt/script
zcat KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
python3 ${script}/summarizeAbundance.py \
      -i KEGG.KO.txt \
        -m ${db}/KEGG/KO1-4.txt \
        -c 2,3,4 -s ',+,+,' -n raw \
        -o KEGG
wc -l KEGG*