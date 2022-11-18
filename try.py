#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 14:36:43 2022

@author: roma
"""
import sys

vcf_file = "/home/roma/Desktop/Korovi/plink/combinded_cow.vcf"
df_file = sys.argv[1]
fam_file = sys.argv[2]
hw_file = sys.argv[3]
ens_file = "/home/roma/Desktop/Ensembl_righr.txt"
import io

import pandas as pd
import numpy as np
import scipy.stats

ens = pd.read_csv(ens_file, sep='\t')
ens = ens[
    ["#Uploaded_variation", "Location", "Consequence", "Gene", "BIOTYPE", "Amino_acids"]
].drop_duplicates()



import pandas as pd 

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
df = pd.read_excel(df_file)

fam = pd.read_csv(fam_file, sep=" ", header=None)
fam = fam[(fam[5] != -9) & (fam[5] != 0)].reset_index(drop=True)
top = 30
df = df.astype({'P': 'float'})



 
df = df.sort_values(by="P").reset_index(drop=True)


hw = pd.read_csv(hw_file, sep="\t")
hw = hw[["SNP", "GENO", "P"]]
hw = hw[hw["SNP"].isin(df["SNP"])].reset_index(drop=True)


final_frame = pd.DataFrame(
    columns=[
        "CHR",
        "SNP",
        "REF",
        "ALT",
        "GT",
        "Consequence",
        "Biotype",
        "Location",
        "Gene",
        "Amino acids",
        "QTL Name",
        "QTL ID",
        "RS",
        "MEAN_PHENO",
        "MEAN CI95",
        "MEDIAN_PHENO",
        "LOWER MEDIAN CI95",
        "UPPER MEDIAN CI95",
        "25th",
        "75th",
        "DELTA",
        "MAF",
        "P Wald test",
        "R-square",
        "Beta",
    ],
    index=range(len(df) * 3),
)
i_v = 0
mv = np.median(np.array(fam[5]))
mean_gs = np.mean(np.array(fam[5]))
for i in df["SNP"]:
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
            if len(fam[fam[0] == t]) != 0:
                v.append(fam[fam[0] == t ][5].to_list()[0])
    for v in v_list:
        final_frame["MAF"][i_v] = df[df["SNP"] == i]["MAF"].to_list()[0]
        final_frame["CHR"][i_v] = df[df["SNP"] == i]["CHR"].to_list()[0]
        final_frame["P Wald test"][i_v] = df[df["SNP"] == i]["P"].to_list()[0]
        final_frame["R-square"][i_v] = df[df["SNP"] == i]["R2"].to_list()[0]
        final_frame["Beta"][i_v] = df[df["SNP"] == i]["BETA"].to_list()[0]
        final_frame["SNP"][i_v] = i
        final_frame["REF"][i_v] = vcf[vcf["ID"] == i].reset_index(drop=True)["REF"][0]
        final_frame["ALT"][i_v] = vcf[vcf["ID"] == i].reset_index(drop=True)["ALT"][0]
        final_frame["Biotype"][i_v] = ",".join(
            map(
                str,
                ens[ens["#Uploaded_variation"] == df[df["SNP"] == i]["RS"].to_list()[0]]["BIOTYPE"]
                .drop_duplicates()
                .to_list(),
            )
        )
        final_frame["Location"][i_v] = ",".join(
            map(
                str,
                ens[ens["#Uploaded_variation"] == df[df["SNP"] == i]["RS"].to_list()[0]]["Location"]
                .drop_duplicates()
                .to_list(),
            )
        )
        final_frame["Amino acids"][i_v] = ",".join(
            map(
                str,
                ens[ens["#Uploaded_variation"] == df[df["SNP"] == i]["RS"].to_list()[0]]["Amino_acids"]
                .drop_duplicates()
                .to_list(),
            )
        )
        final_frame["Gene"][i_v] = ",".join(
            map(
                str,
                ens[ens["#Uploaded_variation"] == df[df["SNP"] == i]["RS"].to_list()[0]]["Gene"]
                .drop_duplicates()
                .to_list(),
            )
        )
        final_frame["Consequence"][i_v] = ",".join(
            map(
                str,
                ens[ens["#Uploaded_variation"] == df[df["SNP"] == i]["RS"].to_list()[0]]["Consequence"]
                .drop_duplicates()
                .to_list(),
            )
        )
        final_frame["QTL ID"][i_v] = df[df["SNP"] == i]["QTL_ID"].to_list()[0]
        final_frame["QTL Name"][i_v] = df[df["SNP"] == i]["QTL_NAME"].to_list()[0]
        final_frame["RS"][i_v] = df[df["SNP"] == i]["RS"].to_list()[0]
        if v == v_w_homo:
            final_frame["GT"][i_v] = (
                vcf[vcf["ID"] == i].reset_index(drop=True)["REF"][0]
                + vcf[vcf["ID"] == i].reset_index(drop=True)["REF"][0]
            )
        elif v == v_hetero:
            final_frame["GT"][i_v] = (
                vcf[vcf["ID"] == i].reset_index(drop=True)["REF"][0]
                + vcf[vcf["ID"] == i].reset_index(drop=True)["ALT"][0]
            )
        else:
            final_frame["GT"][i_v] = (
                vcf[vcf["ID"] == i].reset_index(drop=True)["ALT"][0]
                + vcf[vcf["ID"] == i].reset_index(drop=True)["ALT"][0]
            )
        if len(v) != 0:
            final_frame["MEAN_PHENO"][i_v] = float(mean_confidence_interval(v)[0])
            final_frame["DELTA"][i_v] =  str((float(mean_confidence_interval(v)[0])/mean_gs - 1) * 100) + '%'
            final_frame["LOWER MEDIAN CI95"][i_v] = float(median_confidence_interval(v)[0])
            final_frame["UPPER MEDIAN CI95"][i_v] = float(median_confidence_interval(v)[1])
            final_frame["MEAN CI95"][i_v] = float(mean_confidence_interval(v)[1])
            final_frame["MEDIAN_PHENO"][i_v] = float(mean_confidence_interval(v)[2])
            final_frame["25th"][i_v] = float(mean_confidence_interval(v)[3])
            final_frame["75th"][i_v] = float(mean_confidence_interval(v)[4])
        i_v += 1
final_frame = pd.merge(final_frame, hw, how="left", left_on="SNP", right_on="SNP")
final_frame["TEST"] = ""

for i in range(len(final_frame)):
    if mv < final_frame["LOWER MEDIAN CI95"][i]:
        final_frame["TEST"][i] = "П"
    elif mv > final_frame["UPPER MEDIAN CI95"][i]:
        final_frame["TEST"][i] = "Н"

df_new = pd.DataFrame(columns=final_frame.columns, index = range(1))
final_frame = final_frame.append(df_new).reset_index(drop=True)
print (final_frame.index[-1])
print (final_frame['SNP'][final_frame.index[-1]])
final_frame['SNP'][final_frame.index[-1]] = 'Общее по выборке'
print (final_frame['SNP'][final_frame.index[-1]])
final_frame['MEAN_PHENO'][final_frame.index[-1]] = mean_gs
final_frame['MEDIAN_PHENO'][final_frame.index[-1]] = mv
final_frame['MEAN CI95'][final_frame.index[-1]] = float(mean_confidence_interval(np.array(fam[5]))[1])
final_frame['LOWER MEDIAN CI95'][final_frame.index[-1]] = float(median_confidence_interval(np.array(fam[5]))[0])
final_frame['UPPER MEDIAN CI95'][final_frame.index[-1]] = float(median_confidence_interval(np.array(fam[5]))[1])
final_frame['25th'][final_frame.index[-1]] = float(mean_confidence_interval(np.array(fam[5]))[3])
final_frame['75th'][final_frame.index[-1]] = float(mean_confidence_interval(np.array(fam[5]))[4])
print (final_frame['SNP'][final_frame.index[-1]])



final_frame.to_excel(sys.argv[4], index=False)
