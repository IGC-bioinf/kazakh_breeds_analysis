#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 18:57:51 2022

@author: roma
"""
import sys
file1 = '/home/roma/Desktop/Ensembl_righr.txt'
file2 = sys.argv[1]
file3 = '/home/roma/Desktop/Korovi/plink/AF/af_log.frq'
import pandas as pd
ann_df = pd.read_csv(file1, sep = '\t')
tar_df = pd.read_excel(file2)
frq_df = pd.read_csv(file3, engine='python', sep='\t')


            
df_cd = pd.merge(tar_df, ann_df, how='left', left_on = 'RS', right_on = '#Uploaded_variation')
df_fin = pd.merge(df_cd, frq_df, how='left', left_on = 'SNP', right_on = 'SNP')
#tar_df = df_fin[['SNP', 'BP', 'NMISS', 'BETA', 'SE', 'R2', 'T', 'P', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'MAF', 'EXON', 'Amino_acids']]
tar_df = tar_df.sort_values(by=['SNP']).reset_index(drop=True)
tar_df = tar_df.drop_duplicates(subset=['SNP'], keep = 'last').reset_index(drop=True)
tar_df = tar_df[tar_df['MAF'] >= 0.05]
#tar_df = tar_df[~tar_df['EXON'].fillna('').str.contains('-')].reset_index(drop=True)
#tar_df = tar_df.dropna(subset=['EXON'])   
import numpy as np
            
tar_df['-logp'] = -np.log10(tar_df.P) #tar_df = tar_df.sort_values(['seqnames','start'])
tar_df.reset_index(inplace=True, drop=True)
max_values = tar_df.nlargest(n=10, columns=['BETA'], keep='all')
min_values = tar_df.nsmallest(n=10, columns=['BETA'], keep='all')
max_sum = sum(max_values['BETA'])
min_sum =  sum(min_values['BETA'])
geno = ['0/0', '0/1', '1/1', 'ENS', 'RS', 'REF', 'ALT']
max_tar_df = pd.DataFrame(columns =geno , index = max_values['SNP'])
for i in max_tar_df.index:
    max_tar_df['0/0'][i] = 0
    max_tar_df['0/1'][i] = float(max_values[max_values['SNP'] == i]['BETA'])/max_sum*100
    max_tar_df['1/1'][i] = float(max_values[max_values['SNP'] == i]['BETA']) *2/max_sum*100
    max_tar_df['ENS'][i] = max_values[max_values['SNP'] == i]['Gene'].to_string(index=False)
    max_tar_df['RS'][i] = max_values[max_values['SNP'] == i]['RS'].to_string(index=False)
    max_tar_df['REF'][i] = frq_df[frq_df['SNP'] == i]['A1'].to_string(index=False)
    max_tar_df['ALT'][i] = frq_df[frq_df['SNP'] == i]['A2'].to_string(index=False)
    #max_tar_df['Amino_acids'][i] = max_values[max_values['SNP'] == i]['Amino_acids'].to_string(index=False)
min_tar_df = pd.DataFrame(columns =geno , index = min_values['SNP'])    
for i in min_tar_df.index:
    min_tar_df['0/0'][i] = 0
    min_tar_df['0/1'][i] = float(min_values[min_values['SNP'] == i]['BETA'])/min_sum*100*-1
    min_tar_df['1/1'][i] = float(min_values[min_values['SNP'] == i]['BETA'])  *2/min_sum*100*-1
    min_tar_df['ENS'][i] = min_values[min_values['SNP'] == i]['Gene'].to_string(index=False)
    min_tar_df['RS'][i] = min_values[min_values['SNP'] == i]['RS'].to_string(index=False)
    min_tar_df['REF'][i] = frq_df[frq_df['SNP'] == i]['A1'].to_string(index=False)
    min_tar_df['ALT'][i] = frq_df[frq_df['SNP'] == i]['A2'].to_string(index=False)
    #min_tar_df['Amino_acids'][i] = min_values[min_values['SNP'] == i]['Amino_acids'].to_string(index=False)
    
def write_dataframes_to_excel_sheet(dataframes, dir, name):
    with pd.ExcelWriter(f'{dir}/{name}.xlsx', engine='xlsxwriter') as writer:
        workbook = writer.book
        worksheet = workbook.add_worksheet('Result')
        writer.sheets['Result'] = worksheet

        COLUMN = 0
        row = 0

        for df in dataframes:
            #worksheet.write_string(row, COLUMN)
            #row += 1
            df.to_excel(writer, sheet_name='Result',
                        startrow=row, startcol=COLUMN)
            if len(df.columns) != 1:
                row += df.shape[0] + 2
            else:
                row += df.shape[0] + 3
min_tar_df_zero = min_tar_df.copy()
max_tar_df_zero = max_tar_df.copy()
form = pd.DataFrame(columns=['Значение'], index = [ "Убыль", 'Прибыль', "Итого"])
form['Значение']['Прибыль'] = '=(B2*B14)+(B3*B15)+(B4*B16)+(B5*B17)+(B6*B18)+(B7*B19)+(B8*B20)+(B9*B21)+(B10*B22)+(B11*B23)+(C2*C14)+(C3*C15)+(C4*C16)+(C5*C17)+(C6*C18)+(C7*C19)+(C8*C20)+(C9*C21)+(C10*C22)+(C11*C23)+(D2*D14)+(D3*D15)+(D4*D16)+(D5*D17)+(D6*D18)+(D7*D19)+(D8*D20)+(D9*D21)+(D10*D22)+(D11*D23)'
form['Значение']["Убыль"] = '=(B26*B38)+(B27*B39)+(B28*B40)+(B29*B41)+(B30*B42)+(B31*B43)+(B32*B44)+(B33*B45)+(B34*B46)+(B35*B47)+(C26*C38)+(C27*C39)+(C28*C40)+(C29*C41)+(C30*C42)+(C31*C43)+(C32*C44)+(C33*C45)+(C34*C46)+(C35*C47)+(D26*D38)+(D27*D39)+(D28*D40)+(D29*D41)+(D30*D42)+(D31*D43)+(D32*D44)+(D33*D45)+(D34*D46)+(D35*D47)' 
form['Значение']["Итого"] = '=A49+A50'



#print(min_tar_df.to_csv(), end='\n')
for col in ['0/0', '0/1', '1/1']:
    min_tar_df_zero[col].values[:] = 0
#print(min_tar_df[['0/0', '0/1', '1/1']].to_csv())
#print(max_tar_df.to_csv(), end= '\n')
for col in ['0/0', '0/1', '1/1']:
    max_tar_df_zero[col].values[:] = 0
#print(max_tar_df[['0/0', '0/1', '1/1']].to_csv())


dataframes = [max_tar_df, max_tar_df_zero, min_tar_df, min_tar_df_zero, form]
write_dataframes_to_excel_sheet(dataframes, '/home/roma', sys.argv[2])