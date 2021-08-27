#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:24:53 2021

@author: ryan
"""
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)




with open('RNAcountsBM3D.txt') as f:
    lines = f.readlines()
BM = [[int(j) for j in i.split()] for i in lines[1:]]
tpBM = ['BM' for i in BM]
with open('RNAcountsUC3D.txt') as f:
    lines = f.readlines()
UC = [[int(j) for j in i.split()] for i in lines[1:]]
tpUC = ['UC' for i in UC]
exp = BM+UC
tp=tpBM+tpUC
gene = lines[0].split()

palette = sns.color_palette()

N = 8
E = pd.DataFrame(BM[:N]+UC[-N:],['BM']*N+['UC']*N,gene)
sns.clustermap(E, metric="correlation",
                xticklabels=1,yticklabels=0,
                cmap='icefire',
                # row_cluster=False,col_cluster=False,
                row_colors=[palette[5 + 3*int(i=='UC')] for i in (['BM']*N+['UC']*N)],
                figsize=(7, 5))
plt.savefig('Quantification/TinyExpressionHeatmap3D.png')
E = pd.DataFrame(zscore(BM[:N]+UC[-N:],axis=0),['BM']*N+['UC']*N,gene)
sns.clustermap(E, metric="correlation",
                xticklabels=1,yticklabels=0,
                cmap='icefire',
                # row_cluster=False,col_cluster=False,
                row_colors=[palette[5 + 3*int(i=='UC')] for i in (['BM']*N+['UC']*N)],
                figsize=(7, 5))
plt.savefig('Quantification/TinyZExpressionHeatmap3D.png')