#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:30:48 2022

@author: roma
"""

import os
os.chdir('/home/roma/Desktop/Korovi/figures_for_article')
import matplotlib.ticker as ticker
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt
kbr_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_11_assoc_results.xlsx_withqtl.xlsx'
kbo_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_13_assoc_results.xlsx_withqtl.xlsx'
kbe_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_15_assoc_results.xlsx_withqtl.xlsx'
kbg_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_14_assoc_results.xlsx_withqtl.xlsx'
ar_file = '/home/roma/Desktop/Korovi/final_tables/aul_8_assoc_results.xlsx_withqtl.xlsx'
ao_file = '/home/roma/Desktop/Korovi/final_tables/aul_9_assoc_results.xlsx_withqtl.xlsx'
ae_file = '/home/roma/Desktop/Korovi/final_tables/aul_12_assoc_results.xlsx_withqtl.xlsx'
ag_file = '/home/roma/Desktop/Korovi/final_tables/aul_11_assoc_results.xlsx_withqtl.xlsx'
kbr = pd.read_excel(kbr_file)[pd.read_excel(kbr_file)['CHR'] != 30]
kbo = pd.read_excel(kbo_file)[pd.read_excel(kbo_file)['CHR'] != 30]
kbe = pd.read_excel(kbe_file)[pd.read_excel(kbe_file)['CHR'] != 30]
kbg = pd.read_excel(kbg_file)[pd.read_excel(kbg_file)['CHR'] != 30]
ar = pd.read_excel(ar_file)[pd.read_excel(ar_file)['CHR'] != 30]
ao = pd.read_excel(ao_file)[pd.read_excel(ao_file)['CHR'] != 30]
ae = pd.read_excel(ae_file)[pd.read_excel(ae_file)['CHR'] != 30]
ag = pd.read_excel(ag_file)[pd.read_excel(ag_file)['CHR'] != 30]
min_kbr = 0.05 / len(kbr[kbr['P'] < 0.05])
min_kbo = 0.05 / len(kbo[kbo['P'] < 0.05])
min_kbe = 0.05 / len(kbe[kbe['P'] < 0.05])
min_kbg = 0.05 / len(kbg[kbg['P'] < 0.05])
min_ar = 0.05 / len(ar[ar['P'] < 0.05])
min_ao = 0.05 / len(ao[ao['P'] < 0.05])
min_ae = 0.05 / len(ae[ae['P'] < 0.05])
min_ag = 0.05 / len(ag[ag['P'] < 0.05])
max_kbr = 0.05 / len(kbr)
max_kbo = 0.05 / len(kbo)
max_kbe = 0.05 / len(kbe)
max_kbg = 0.05 / len(kbg)
max_ar = 0.05 / len(ar)
max_ao = 0.05 / len(ao)
max_ae = 0.05 / len(ae)
max_ag = 0.05 / len(ag)
values_list = [min_kbr, min_kbo, min_kbe, min_kbg, min_ar, min_ao, min_ae, min_ag, max_kbr, max_kbo, max_kbe, max_kbg, max_ar, max_ao, max_ae, max_ag]



kbr[kbr['P'] < min_kbr].to_excel(kbr_file + '_min.xlsx', index = False)
kbo[kbo['P'] < min_kbo].to_excel(kbo_file + '_min.xlsx', index = False)
kbe[kbe['P'] < min_kbe].to_excel(kbe_file + '_min.xlsx', index = False)
kbg[kbg['P'] < min_kbg].to_excel(kbg_file + '_min.xlsx', index = False)
ar[ar['P'] < min_ar].to_excel(ar_file + '_min.xlsx', index = False)
ao[ao['P'] < min_ao].to_excel(ao_file + '_min.xlsx', index = False)
ae[ae['P'] < min_ae].to_excel(ae_file + '_min.xlsx', index = False)
ag[ag['P'] < min_ag].to_excel(ag_file + '_min.xlsx', index = False)


kbr[kbr['P'] < max_kbr].to_excel(kbr_file + '_max.xlsx', index = False)
kbo[kbo['P'] < max_kbo].to_excel(kbo_file + '_max.xlsx', index = False)
kbe[kbe['P'] < max_kbe].to_excel(kbe_file + '_max.xlsx', index = False)
kbg[kbg['P'] < max_kbg].to_excel(kbg_file + '_max.xlsx', index = False)
ar[ar['P'] < max_ar].to_excel(ar_file + '_max.xlsx', index = False)
ao[ao['P'] < max_ao].to_excel(ao_file + '_max.xlsx', index = False)
ae[ae['P'] < max_ae].to_excel(ae_file + '_max.xlsx', index = False)
ag[ag['P'] < max_ag].to_excel(ag_file + '_max.xlsx', index = False)







a_list = [kbr, kbo , kbe, kbg,  ar, ao, ae, ag]
n_list = ['Kazakh white-headed Birth weight', 'Kazakh white-headed Weaning weight' , 'Kazakh white-headed Average daily gain', 'Kazakh white-headed Yearling weight', 'Auliekol Birth weight', 'Auliekol Weaning weight', 'Auliekol Average daily gain', 'Auliekol Yearling weight']
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import pylab as py
for i, n in zip(a_list, n_list):
    min_v = 0.05 / len(i[i['P'] < 0.05])
    max_v = 0.05/len(i)
    i['-logp'] = -np.log10(i.P); i = i.sort_values(['CHR','BP'])
    i.reset_index(inplace=True, drop=True); i['i'] = i.index
    
    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    plot = sns.relplot(data=i, x='i', y='-logp', aspect=3.7, 
                       hue='CHR', palette = 'bright', legend=None) 
    chrom_i=i.groupby('CHR')['i'].median()
    plot.ax.set_xlabel('Chromosome'); plot.ax.set_xticks(chrom_i);
    plot.ax.set_xticklabels(chrom_i.index)
    plot.fig.suptitle(n)
    plt.plot(i['i'], [-np.log10(min_v)] * len(i['i']))
    plt.plot(i['i'], [-np.log10(max_v)] * len(i['i']))
    plt.ylim(0,9)
    plot.savefig('Manhattan plot_' + n +'.png',dpi=100,  bbox_inches = 'tight') ###################
    print('PLOT DONE')
   # sm.qqplot(i['P'], line ='45')
    #py.show()


thresh = 0.01118
kbr = kbr[kbr['P'] < thresh]
kbo = kbo[kbo['P'] < thresh]
kbe = kbe[kbe['P'] < thresh]
kbg = kbg[kbg['P'] < thresh]
ar = ar[ar['P'] < thresh]
ao = ao[ao['P'] < thresh]
ae = ae[ae['P'] < thresh]
ag = ag[ag['P'] < thresh]
k_list = [set(kbr['SNP'].to_list()), set(kbo['SNP'].to_list()), set(kbg['SNP'].to_list()),set(ar['SNP'].to_list()), set(ao['SNP'].to_list()), set(ag['SNP'].to_list())]
a_list = [kbr, kbo , kbe, kbg, ar, ao, ae, ag]
df = pd.DataFrame(columns=['Kazakh white-headed', 'Auliekol', 'Kazakh white-headed common', 'Auliekol common'], index=['Birth weight', 'Weaning weight', "Yearling weight", 'Average daily gain'])
df['Kazakh white-headed']['Birth weight'] = len(kbr)
df['Kazakh white-headed']['Weaning weight'] = len(kbo)
df['Kazakh white-headed']['Average daily gain'] = len(kbe)
df['Kazakh white-headed']['Yearling weight'] = len(kbg)
df['Auliekol']['Birth weight'] = len(ar)*-1
df['Auliekol']['Weaning weight'] = len(ao)*-1
df['Auliekol']['Average daily gain'] = len(ae)*-1
df['Auliekol']['Yearling weight'] = len(ag)*-1
df['Kazakh white-headed common']['Birth weight'] = len(list(set(kbr.SNP) & set(ar.SNP)))
df['Kazakh white-headed common']['Weaning weight'] = len(list(set(kbo.SNP) & set(ao.SNP)))
df['Kazakh white-headed common']['Average daily gain'] = len(list(set(kbe.SNP) & set(ae.SNP)))
df['Kazakh white-headed common']['Yearling weight'] = len(list(set(kbg.SNP) & set(ag.SNP)))
df['Auliekol common']['Birth weight'] = df['Kazakh white-headed common']['Birth weight']*-1
df['Auliekol common']['Weaning weight'] = df['Kazakh white-headed common']['Weaning weight']*-1
df['Auliekol common']['Average daily gain'] = df['Kazakh white-headed common']['Average daily gain']*-1
df['Auliekol common']['Yearling weight'] = df['Kazakh white-headed common']['Yearling weight']*-1
   

def neg_to_pos(n, position):
    n = int(n)
    if n < 0:
        return str(n * -1)
    else:
        return str(n)

max_v = df.abs().max().max()
fig, ax = plt.subplots(1, figsize=(18, 5))

plt.barh(df.index, df[df.columns[0]], color = '#337AE3', height = 0.5, align='edge')
plt.barh(df.index, df[df.columns[1]], color = '#DB4444', height = 0.5, align='edge')
plt.barh(df.index, df[df.columns[3]], color = 'green', height = 0.5, align='edge')
plt.barh(df.index, df[df.columns[2]], color = 'green', height = 0.5, align='edge')


plt.xlim((max_v*-1)-10, max_v+10)

ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_axisbelow(True)
#px= ax.twiny()
ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
#ax.xaxis.set_major_formatter(ticker.PercentFormatter())
ax.xaxis.set_major_formatter(ticker.FuncFormatter(neg_to_pos))

   
plt.xticks(rotation=90)
ax.xaxis.grid(color='gray', linestyle='dashed', alpha=0.7)
ax.tick_params(axis='both', which='major', labelsize=15)
#plt.xticks(df.index)
legend_label = list(df.columns)
plt.legend(legend_label, ncol=2, frameon = False, prop={'size': 15}, bbox_to_anchor=(1.04,1), loc="upper left")

plt.savefig('minus0.01_bar.png', dpi=100,  bbox_inches = 'tight')
plt.show()


os.chdir('/home/roma/Desktop/Korovi')
from pyvenn import venn
os.chdir('/home/roma/Desktop/Korovi/figures_for_article')
labels = venn.get_labels(k_list, fill=['number'])
fig, ax = venn.venn6(labels, names=['Kazakh white-headed Birth weight', 'Kazakh white-headed Weaning weight', 'Kazakh white-headed Yearling weight', 'Auliekol Birth weight', 'Auliekol Weaning weight',  'Auliekol Yearling weight'])
fig.show()
plt.savefig('minus0.01_venn.png', dpi=100,  bbox_inches = 'tight')



chr_df = pd.DataFrame(columns=['Kazakh white-headed Birth weight', 'Kazakh white-headed Weaning weight', 'Kazakh white-headed Average daily gain', 'Kazakh white-headed Yearling weight', 'Auliekol Birth weight', 'Auliekol Weaning weight', 'Auliekol Average daily gain', 'Auliekol Yearling weight'], index = range(1,31) )
for a,n in zip(a_list, n_list):
   for k in chr_df.index:
       chr_df[n][k] = float(len(a[a['CHR'] == k]))


chr_df = chr_df.astype(float)
chr_df1 = chr_df[['Kazakh white-headed Birth weight', 'Kazakh white-headed Weaning weight',  'Kazakh white-headed Yearling weight', 'Kazakh white-headed Average daily gain']]
chr_df2 = chr_df[['Auliekol Birth weight', 'Auliekol Weaning weight', 'Auliekol Yearling weight', 'Auliekol Average daily gain' ]]

sns.set(rc = {'figure.figsize':(15,8)})
svm1 = sns.heatmap(chr_df1, cmap="YlGnBu")
figure1 = svm1.get_figure()    
figure1.savefig('kaz_bel_heat_minus0.01.png', dpi=400,  bbox_inches = 'tight')
sns.set(rc = {'figure.figsize':(15,8)})
svm2 = sns.heatmap(chr_df2, cmap="YlGnBu")
figure2 = svm2.get_figure()    
figure2.savefig('aul_heat_minus0.01.png', dpi=400,  bbox_inches = 'tight')


#from qmplot import qqplot

#ax = qqplot(data=ar["P"], figname="output_QQ_plot.png", dpi=300)







