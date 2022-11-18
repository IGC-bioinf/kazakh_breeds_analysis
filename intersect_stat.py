#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 04:33:12 2022

@author: roma
"""
import os
os.chdir('/home/roma/Desktop/Korovi/figures_for_article')
import matplotlib.ticker as ticker
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt
kbr_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_11_add_stat.xlsx'
kbo_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_13_add_stat.xlsx'
kbe_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_15_add_stat.xlsx'
kbg_file = '/home/roma/Desktop/Korovi/final_tables/kaz_bel_14_add_stat.xlsx'
ar_file = '/home/roma/Desktop/Korovi/final_tables/aul_8_add_stat.xlsx'
ao_file = '/home/roma/Desktop/Korovi/final_tables/aul_9_add_stat.xlsx'
ae_file = '/home/roma/Desktop/Korovi/final_tables/aul_12_add_stat.xlsx'
ag_file = '/home/roma/Desktop/Korovi/final_tables/aul_11_add_stat.xlsx'
kbr = pd.read_excel(kbr_file)
kbo = pd.read_excel(kbo_file)
kbe = pd.read_excel(kbe_file)
kbg = pd.read_excel(kbg_file)
ar = pd.read_excel(ar_file)
ao = pd.read_excel(ao_file)
ae = pd.read_excel(ae_file)
ag = pd.read_excel(ag_file)

keys_frames = {"КБ ЖМР" : kbr, "КБ ЖМО": kbo, "КБ ЕП" : kbe, "КБ ЖМГ" : kbg, "АУЛ ЖМР" : ar, "АУЛ ЖМО": ao, "АУЛ ЕП" : ae, "АУЛ ЖМГ" : ag}

for i in list(keys_frames.keys()):
    tmp_list = [s for s in list(keys_frames.keys()) if s != i]
    keys_frames[i]['OTHER TOP30'] = ''
    for k in range(len(keys_frames[i])):
        for l in tmp_list:
            if len(keys_frames[l][keys_frames[l]['SNP'] == keys_frames[i]['SNP'][k]]) != 0:
                keys_frames[i]['OTHER TOP30'][k] = keys_frames[i]['OTHER TOP30'][k] + '|' + l
                

tkbr_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/kaz_bel_11_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tkbo_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/kaz_bel_13_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tkbe_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/kaz_bel_15_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tkbg_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/kaz_bel_14_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tar_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/aul_8_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tao_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/aul_9_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tae_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/aul_12_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'
tag_file = '/home/roma/Desktop/Korovi/final_tables/Shablons/aul_11_assoc_results.xlsx_withqtl_add_inf.xlsx.xlsx'                
tkbr = pd.read_excel(tkbr_file)
tkbo = pd.read_excel(tkbo_file)
tkbe = pd.read_excel(tkbe_file)
tkbg = pd.read_excel(tkbg_file)
tar = pd.read_excel(tar_file)
tao = pd.read_excel(tao_file)
tae = pd.read_excel(tae_file)
tag = pd.read_excel(tag_file)
tkeys_frames = {"КБ ЖМР" : tkbr, "КБ ЖМО": tkbo, "КБ ЕП" : tkbe, "КБ ЖМГ" : tkbg, "АУЛ ЖМР" : tar, "АУЛ ЖМО": tao, "АУЛ ЕП" : tae, "АУЛ ЖМГ" : tag}
for i in list(tkeys_frames.keys()):
    keys_frames[i]['OTHER TOP10'] = ''
    for k in range(len(keys_frames[i])):
        for l in list(tkeys_frames.keys()):
            if len(tkeys_frames[l][tkeys_frames[l]['SNP'] == keys_frames[i]['SNP'][k]]) != 0:
                keys_frames[i]['OTHER TOP10'][k] = keys_frames[i]['OTHER TOP10'][k] + '|' + l
kbr.to_excel(kbr_file+'_stat.xlsx', index = False)
kbo.to_excel(kbo_file+'_stat.xlsx', index = False)
kbe.to_excel(kbe_file+'_stat.xlsx', index = False)
kbg.to_excel(kbg_file+'_stat.xlsx', index = False)
ar.to_excel(ar_file+'_stat.xlsx', index = False)
ao.to_excel(ao_file+'_stat.xlsx', index = False)
ae.to_excel(ae_file+'_stat.xlsx', index = False)
ag.to_excel(ag_file+'_stat.xlsx', index = False)
                
                

                
                
