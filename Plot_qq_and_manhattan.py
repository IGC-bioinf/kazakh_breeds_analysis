#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 14:05:46 2022

@author: roma
"""
from qmplot import qqplot
from qmplot import manhattanplot
import sys
import pandas as pd
file = '/home/roma/Desktop/Korovi/final_tables/aul_12_assoc_results.xlsx_withqtl.xlsx'
df = pd.read_excel(file)

ax = qqplot(data=df["P"], figname="output_qq_plot.png", dpi=300)
ax = manhattanplot(data=df, figname="output_manhattan_plot.png")

import numpy as np
import seaborn as sns
df['-logp'] = -np.log10(df.P); df = df.sort_values(['CHR','BP'])
df.reset_index(inplace=True, drop=True); df['i'] = df.index

# Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
plot = sns.relplot(data=df, x='i', y='-logp', aspect=3.7, 
                   hue='CHR', palette = 'bright', legend=None) 
chrom_df=df.groupby('CHR')['i'].median()
plot.ax.set_xlabel('Chromosome'); plot.ax.set_xticks(chrom_df);
plot.ax.set_xticklabels(chrom_df.index)
#plot.fig.suptitle(sys.argv[3])
plot.savefig('Manhattan plot_' + sys.argv[3] +'.png',dpi=100,  bbox_inches = 'tight') ###################
print('PLOT DONE')