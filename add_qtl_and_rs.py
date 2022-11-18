#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 17:14:46 2022

@author: roma
"""
import pandas as pd
import sys
gff_file = '/home/roma/Desktop/Korovi/qtl/QTLdb_cattleUMD3.gff'
file1 = '/home/roma/Desktop/Korovi/qtl/rs_id.txt'
file2 = sys.argv[1]
df_rs = pd.read_csv(file1, names=['RS', 'SNP'])
df_ann = pd.read_excel(file2)


df_f = pd.merge(df_ann, df_rs, how='left', left_on = 'SNP', right_on = 'SNP')
def readGff(path):
    with open(path,'r', errors='replace') as gff_file: 
        gff_lines = gff_file.readlines()
    gff_data = [line.replace('\t', '=').replace(';', '=').split('=') for line in gff_lines if '#' != line[0]]
    gff_frame = pd.DataFrame(gff_data)
    gff_frame = gff_frame[[9,11,22]]
    gff_frame = gff_frame.rename(columns={9: "QTL_ID", 11: "QTL_NAME", 22: 'RS'})
    return gff_frame
gff = readGff(gff_file)
gff['RS']=gff['RS'].str.split(",")
gff = gff.explode('RS').reset_index(drop=True)
df_f = pd.merge(df_f, gff, how='left', left_on = 'RS', right_on = 'RS')
ensenmbl_file = '/home/roma/Desktop/Ensembl_righr.txt'
df_ens =  pd.read_csv(ensenmbl_file, sep='\t')
df_f = pd.merge(df_f, df_ens[['#Uploaded_variation', 'Location', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene']], how='left', left_on = 'RS', right_on = '#Uploaded_variation').drop_duplicates().reset_index(drop=True)

file3 = '/home/roma/Desktop/Korovi/plink/AF/af_log.frq'
frq_df = pd.read_csv(file3, engine='python', sep='\t')
df_f = pd.merge(df_f, frq_df[['SNP', 'MAF']], how='left', left_on = 'SNP', right_on = 'SNP')
del df_f['#Uploaded_variation']
for i in df_f.columns:    
    df_f[i] = df_f[i].apply(str)
df_f = df_f.groupby(['SNP', 'Location', 'BETA', 'P', 'CHR', 'BP','MAF', 'T', 'R2', 'SE', 'RS']).agg({'QTL_ID': ','.join, 'QTL_NAME': ','.join, 'Consequence': ','.join, 'IMPACT': ','.join, 'SYMBOL': ','.join, 'Gene': ','.join}) 
for i in ['QTL_ID', 'QTL_NAME', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene']:
    for k in range(len(df_f[i])):
        df_f[i][k] = str(list(set(df_f[i][k].split(',')))).replace('[', '').replace(']', '').replace("'", '')
print('Saving to file ' + sys.argv[2])
df_f.to_excel(sys.argv[2], index=True)