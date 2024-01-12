import pandas as pd
df = pd.read_csv('/workdir/zhangml/20221024_NTM/restart/3.pfam/2.matrix/matrix_pfam_335_growth_numeric.csv', index_col=0)
seq_totals = df.sum(axis=0)
proportions = df.divide(seq_totals, axis=1)
proportions.to_csv('/workdir/zhangml/20221024_NTM/restart/14.cog_kegg_heatmap_new/pfam/gene_proportions.csv')
