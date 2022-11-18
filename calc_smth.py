#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 14:47:43 2022

@author: roma
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
file = sys.argv[1]
df = pd.read_excel(file)
df = df.sort_values(by=['BETA','P']).reset_index(drop=True)
#plt.plot(df['BETA'], df['P'])
df['-logp'] = -np.log10(df.P) #df = df.sort_values(['seqnames','start'])
df.reset_index(inplace=True, drop=True)
max_values = df.nlargest(n=10, columns=['BETA'], keep='all')
min_values = df.nsmallest(n=10, columns=['BETA'], keep='all')
max_sum = sum(max_values['BETA'])
min_sum =  sum(min_values['BETA'])
geno = ['0/0', '0/1', '1/1']
max_df = pd.DataFrame(columns =geno , index = max_values['SNP'])
for i in max_df.index:
    max_df['0/0'][i] = 0
    max_df['0/1'][i] = float(max_values[max_values['SNP'] == i]['BETA'])/max_sum*100
    max_df['1/1'][i] = float(max_values[max_values['SNP'] == i]['BETA']) *2/max_sum*100
min_df = pd.DataFrame(columns =geno , index = min_values['SNP'])    
for i in min_df.index:
    min_df['0/0'][i] = 0
    min_df['0/1'][i] = float(min_values[min_values['SNP'] == i]['BETA'])/min_sum*100*-1
    min_df['1/1'][i] = float(min_values[min_values['SNP'] == i]['BETA']) *2/min_sum*100*-1

print(min_df.to_csv(), end='\n')
for col in min_df.columns:
    min_df[col].values[:] = 0
print(min_df.to_csv())
print(max_df.to_csv(), end= '\n')
for col in max_df.columns:
    max_df[col].values[:] = 0
print(max_df.to_csv())

k1=['B','C','D']
k2 = ['G','H','I']
n1 = range(2,12)
n2 = range(26,36)
k1_str =''
k2_str = ''

for k,z in zip(k1, k2):
    for n in n1:
        k1_str += '(' + k + str(n) + '*' + k+ str(n+12)+ ')' + '+'
        

for k,z in zip(k1, k2):
    for n in n2:
        k2_str += '(' + k + str(n) + '*' + k+ str(n+12)+ ')' + '+'