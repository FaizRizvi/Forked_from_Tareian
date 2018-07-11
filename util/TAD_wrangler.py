import pandas as pd
import numpy as np
import seaborn as sns

cd /Users/caz3so/Dropbox/thesis/data/qATAC_pipeline/tad_gene/

#column types

header_names =['Chromsome_A', 'Start_A', 'Stop_A', 'Gene', 'NA1' , 'Strand', 'TSS_start', 'TSS_stop', 'NA2', 'NA3', 'Exon1', 'Exon2', 'Chromsome_B', 'Start_B', 'Stop_B', 'TAD_name']

# Create a dataframe with the four feature variables
df = pd.read_csv('tad_gene_1bp.bed', header = None, sep='\t', names=header_names)

subset = df[['TAD_name', 'Gene']]

ax = sns.countplot(x='TAD_name', data=df, order=df['TAD_name'].value_counts().index) 

                 data=df_faculty, color="grey", order=df_faculty['PhdDecade'].value_counts().index)
                 
subset = df.to_csv('TAD_Gene_list.txt')