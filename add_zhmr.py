#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 18:47:39 2022

@author: roma
"""
ep_file = '/home/roma/Desktop/Korovi/final_tables/aul_12_add_stat.xlsx_stat.xlsx'
import pandas as pd

import io
import os
import pandas as pd
import numpy as np
import scipy.stats



ep = pd.read_excel(ep_file)
famr_file = '/home/roma/Desktop/Korovi/final_tables/aul_8.fam'
famg_file = '/home/roma/Desktop/Korovi/final_tables/aul_11.fam'

famr = pd.read_csv(famr_file, sep = ' ', header=None)
famr = famr[(famr[5] != -9) & (famr[5] != 0)]
famg = pd.read_csv(famg_file, sep = ' ', header=None)
famg = famg[(famg[5] != -9) & (famg[5] != 0)]
ep['MEAN ЖМР'] = ""
ep['MEAN ЖМГ'] = ""


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m = np.mean(a)
    median, se = np.median(a), scipy.stats.sem(a)
    first = np.percentile(a, 25)
    third = np.percentile(a, 75)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2.0, n - 1)
    return m, h, median, first, third
import io
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

vcf_file = "/home/roma/Desktop/Korovi/plink/combinded_cow.vcf"
vcf = read_vcf(vcf_file)

for i in list(set(ep['SNP'].to_list()[:-1])):
    if ('Общее по выборке' not in str(i)):
        if str(i) != 'nan':
            s_w_homo, s_hetero, s_m_homo = [], [], []
        
            tmp_frame = vcf[vcf["ID"] == i].reset_index(drop=True)
            for k in tmp_frame.columns:
                if "0/0" in str(tmp_frame[k][0]):
                    s_w_homo.append(str(k))
                elif "0/1" in str(tmp_frame[k][0]):
                    s_hetero.append(str(k))
                elif "1/1" in str(tmp_frame[k][0]):
                    s_m_homo.append(str(k))
            s_list = [s_w_homo, s_hetero, s_m_homo]
            v_w_homo, v_hetero, v_m_homo = [], [], []
            v_list = [v_w_homo, v_hetero, v_m_homo]
            
            for s, v in zip(s_list, v_list):
                for t in s:
                    if len(famg[famg[0] == t]) != 0:
                        v.append(famg[famg[0] == t ][5].to_list()[0])
            i_v = 0 
            for v in v_list:
                if len(v) != 0:
                    ep["MEAN ЖМГ"][ep[ep['SNP'] == i].index[0+i_v]] = float(mean_confidence_interval(v)[0])
                    i_v +=1
                
            print (len(v_list))
        
        
    

for i in list(set(ep['SNP'].to_list()[:-1])):
    if ('Общее по выборке' not in str(i)):
        if str(i) != 'nan':
            s_w_homo, s_hetero, s_m_homo = [], [], []
        
            tmp_frame = vcf[vcf["ID"] == i].reset_index(drop=True)
            for k in tmp_frame.columns:
                if "0/0" in str(tmp_frame[k][0]):
                    s_w_homo.append(str(k))
                elif "0/1" in str(tmp_frame[k][0]):
                    s_hetero.append(str(k))
                elif "1/1" in str(tmp_frame[k][0]):
                    s_m_homo.append(str(k))
            s_list = [s_w_homo, s_hetero, s_m_homo]
            v_w_homo, v_hetero, v_m_homo = [], [], []
            v_list = [v_w_homo, v_hetero, v_m_homo]
            for s, v in zip(s_list, v_list):
                for t in s:
                    if len(famr[famr[0] == t]) != 0:
                        v.append(famr[famr[0] == t ][5].to_list()[0])
    
            i_v = 0 
            for v in v_list:
                if len(v) != 0:
                   ep["MEAN ЖМР"][ep[ep['SNP'] == i].index[0+i_v]] = float(mean_confidence_interval(v)[0])
                   i_v +=1
                   
            print (len(v_list))
            
ep.to_excel('/home/roma/Desktop/Korovi/final_tables/test.xlsx', index=False)