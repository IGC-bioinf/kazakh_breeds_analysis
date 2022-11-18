#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 20:49:48 2022

@author: roma
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy.stats
import io
import pandas as pd
tmp_frame=pd.DataFrame()

vcf_file = '/home/roma/Desktop/Korovi/final_tables/stat_final/to_pair/SNP_kaz_bel_zhmg.list.recode.vcf'
fam_file = '/home/roma/Desktop/Korovi/final_tables/stat_final/to_pair/kaz_bel_14.fam'
info_file = '/home/roma/Desktop/Korovi/final_tables/stat_final/to_pair/kaz_bel_14_add_stat.xlsx_stat.xlsx'
savepath = '/home/roma/Desktop/Korovi/new_box_plots/kaz_bel_14/'

file1 = '/home/roma/Desktop/Ensembl_righr.txt'
ann_df = pd.read_csv(file1, sep = '\t')


def median_confidence_interval(data, ci=0.95, p=0.5):
    data = 1.0 * np.array(data)
    if type(data) is pd.Series or type(data) is pd.DataFrame:
        data = data.values

	#flat to one dimension array
    data = data.reshape(-1)
    data = np.sort(data)
    N = data.shape[0]
    lowCount, upCount = scipy.stats.binom.interval(ci, N, p, loc=0)
    lowCount -= 1
    upCount -= 1
    if lowCount < 0:
        lowCount = 0
    if upCount > len(data):
        upCount = len(data)
    
	# print lowCount, upCount
    return data[int(lowCount)], data[int(upCount)]




def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m = np.mean(a)
    median, se = np.median(a), scipy.stats.sem(a)
    first = np.percentile(a, 25)
    third = np.percentile(a, 75)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2.0, n - 1)
    return m, h, median, first, third


def read_vcf(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    return pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})

vcf = read_vcf(vcf_file)
list_of_SNP = vcf['ID'].to_list()
for i in vcf.columns:
    for k in range(len(vcf)):
        if vcf[i][k] == '0/0':
            vcf[i][k] = str(vcf['REF'][k]+'/'+vcf['REF'][k])
        elif vcf[i][k] == '0/1':
            vcf[i][k] = str(vcf['REF'][k]+'/'+vcf['ALT'][k])
        elif vcf[i][k] == '1/1':
            vcf[i][k] = str(vcf['ALT'][k]+'/'+vcf['ALT'][k])
            
result = []
import itertools
pairs = itertools.permutations(list_of_SNP,2)
for i in pairs:
    result.append(i)

nr = list(set((a,b) if a<=b else (b,a) for a,b in result))


fam  = pd.read_csv(fam_file, sep = ' ') 

info = pd.read_excel(info_file)

figure, axis = plt.subplots(15, 3, figsize=(40,60))  
x=0
y=0
tit_num = 1
for pair in nr:
    tmp_vcf = vcf[vcf['ID'].str.contains(pair[0]) | vcf['ID'].str.contains(pair[1]) ]
    fam =  pd.read_csv(fam_file, sep = ' ', header=None)
    columns =  fam[0].to_list()
    to_app = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    for i in to_app:  
        columns.append(i)
    tmp_vcf = tmp_vcf[columns]
    geno_list = []
    for i in tmp_vcf[fam[0]].columns:
        geno_list.append(tmp_vcf[i].to_list())
    
    unique_data = [list(x) for x in set(tuple(x) for x in geno_list)]
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    for i in unique_data:
            if '.' in str(i) or './.' in str(i):
                unique_data.remove(i)
    
    
    values_frame = pd.DataFrame(columns = [str(x) for x in unique_data], index =range(len(columns)))
    
    i_v=0
    
    for i in unique_data:
        for k in tmp_vcf[fam[0]].columns:
            if tmp_vcf[k].to_list() == i:
                if float(fam[fam[0] == k][5]) != -9 and float(fam[fam[0] == k][5]) != 0:    
                    values_frame[str(i)][i_v] = float(fam[fam[0] == k][5])
                    i_v+=1
        i_v = 0
        
    
    values_frame['ALL'] = ''
    for i in range(len(fam)):
        if float(fam[5][i]) != -9 and float(fam[5][i]) != 0:
            values_frame['ALL'][i] = float(fam[5][i])
    
    
    values_frame = values_frame.fillna('')
    
    values_dict = {}
    for i in values_frame.columns:
        values_dict[i] = []
        
    for i in values_frame.columns:
        for k in values_frame[i]:
            tmp = values_dict[i]
            if k != '':
                tmp.append(k)
                values_dict[i] = tmp
                
    import collections
    values_dict = collections.OrderedDict(sorted(values_dict.items()))
    stat_frame = pd.DataFrame(columns = ['ID', 'Mean', 'Median', 'UCI', 'LCI', '25th', '75th', 'Delta', 'Freq'], index = list(values_dict.keys()))
    for i in stat_frame.index:
       
            if i != 'ALL':
                stat_frame['ID'][i] = str(pair[1]) + ' + ' + str(pair[0])
            if len(values_dict[i]) != 0 :
                stat_frame['Mean'][i] = mean_confidence_interval(values_dict[i])[0]
                stat_frame['Median'][i] = mean_confidence_interval(values_dict[i])[2]
                stat_frame['UCI'][i] = median_confidence_interval(values_dict[i])[1]
                stat_frame['LCI'][i] = median_confidence_interval(values_dict[i])[0]
                stat_frame['25th'][i] =mean_confidence_interval(values_dict[i])[3]
                stat_frame['75th'][i] =mean_confidence_interval(values_dict[i])[4]
                stat_frame['Freq'][i] = len(values_dict[i])
    for i in stat_frame.index:
        if i != 'ALL':
            stat_frame['Delta'][i] = stat_frame['Mean'][i] - stat_frame['Mean']['ALL']
        else:
            stat_frame['Delta'][i] = 0
    
    
    to_del=[]
    for k in values_dict.keys():
        if  (stat_frame['Median'][k] < stat_frame['UCI']["ALL"] or  stat_frame['Freq'][k] < 10) and k != 'ALL':
            to_del.append(k)
    for k in to_del:
        del values_dict[k]
    # Python 3.5+
    labels, data = [*zip(*values_dict.items())]  # 'transpose' items to parallel key, value lists
    
    # or backwards compatable    
    labels, data = values_dict.keys(), values_dict.values()
    #Creating subplot of each column with its own scale
    #fig, ax = plt.subplots(1, figsize=(15, 10))
    boxes = axis[x, y].boxplot(data, patch_artist=True)
    for box in boxes["boxes"]:
        box.set(facecolor = "orange")
    axis[x, y].set_xticklabels( labels, rotation=0)
    axis[x, y].set_title(tit_num, loc='center')
    #red_circle = dict(markerfacecolor='red', marker='o', markeredgecolor='white')
    #flierprops = dict(marker='*', markerfacecolor='green', markersize=12,
    #              linestyle='none')
    #fig, axs = plt.subplots(1, len(labels), figsize=(20,10))
    #fig.ylim(0,400)
    #for i, ax in enumerate(axs.flat):
    #        ax.set_ylim(0,400)
    #        ax.boxplot(values_dict[list(labels)[i]], flierprops=red_circle)
    #        ax.set_title(list(values_dict.keys())[i], fontsize=10, fontweight='bold')
    #        ax.tick_params(axis='y', labelsize=12)

    #axis[x,y].tight_layout()
    stat_frame['1st ENS'] = ''
    stat_frame['1st RS'] = ''
    stat_frame['2st ENS'] = ''
    stat_frame['2st RS'] = ''
    stat_frame['1st Symbol'] = ''
    stat_frame['2st Symbol'] = ''
    for ind in stat_frame.index:
        if ind != 'ALL':
            stat_frame['1st ENS'][ind] = info[info['SNP'] == stat_frame['ID'][ind].split('+')[0].replace(' ', '')]['Gene'].drop_duplicates().to_list()[0]
            stat_frame['1st RS'][ind] = info[info['SNP'] == stat_frame['ID'][ind].split('+')[0].replace(' ', '')]['RS'].drop_duplicates().to_list()[0]
            stat_frame['2st ENS'][ind] = info[info['SNP'] == stat_frame['ID'][ind].split('+')[1].replace(' ', '')]['Gene'].drop_duplicates().to_list()[0]
            stat_frame['2st RS'][ind] = info[info['SNP'] == stat_frame['ID'][ind].split('+')[1].replace(' ', '')]['RS'].drop_duplicates().to_list()[0]
            stat_frame['1st Symbol'][ind] = ann_df[ann_df['Gene'] == stat_frame['1st ENS'][ind]]['SYMBOL'].drop_duplicates().to_list()[0]
            stat_frame['2st Symbol'][ind] = ann_df[ann_df['Gene'] == stat_frame['2st ENS'][ind]]['SYMBOL'].drop_duplicates().to_list()[0]
    stat_frame.to_excel(savepath +str(tit_num) + '___'+ str(pair[1]) + ' + ' + str(pair[0]) + '.xlsx')
    tit_num +=1
    if x!=14:
        x+=1
    else:
        x=0
        y+=1
    if len(tmp_frame) != 0:
        tmp_frame = pd.concat([tmp_frame, stat_frame])
    else:
        tmp_frame = stat_frame.copy()
    
#plt.savefig('./box_plots/' + str(pair[0]) + ' + ' + str(pair[1]) + '.png', dpi=300,  bbox_inches = 'tight') 
#plt.subplots_adjust (hspace = 5)
plt.savefig(savepath + 'plot.png', dpi=300,  bbox_inches = 'tight') 
plt.show()
tmp_frame = tmp_frame.drop_duplicates()
tmp_frame.to_excel(savepath + 'comdined.xlsx')   
    

            
            