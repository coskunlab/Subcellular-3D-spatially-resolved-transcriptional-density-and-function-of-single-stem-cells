#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:25:26 2021

@author: ryan
"""
import seaborn as sns
import matplotlib
import pandas as pd
from scipy.stats import zscore
import matplotlib.pyplot as plt
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)




for samp in ['','3D']:
    with open('RNAcountsBM'+samp+'.txt') as f:
        lines = f.readlines()
    BM = [[int(j) for j in i.split()] for i in lines[1:]]
    tpBM = ['BM' for i in BM]
    with open('RNAcountsUC'+samp+'.txt') as f:
        lines = f.readlines()
    UC = [[int(j) for j in i.split()] for i in lines[1:]]
    tpUC = ['UC' for i in UC]
    exp = BM+UC
    tp=tpBM+tpUC
    gene = lines[0].split()
    
    palette = sns.color_palette()
    
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 32}
    matplotlib.rc('font', **font)
    E = pd.DataFrame(exp,tp,gene)
    sns.clustermap(E, metric="correlation",
                    xticklabels=1,yticklabels=0,
                    cmap='icefire',
                    # dendrogram_ratio=(0.2,0.001),
                    # row_cluster=False,col_cluster=False,
                    row_colors=[palette[5 + 3*int(i=='UC')] for i in tp],
                    figsize=(7, 70))
    plt.savefig('Quantification/ExpressionHeatmap'+samp+'.png',dpi=936)
    E = pd.DataFrame(zscore(exp,axis=0),tp,gene)
    sns.clustermap(E, metric="correlation",
                    xticklabels=1,yticklabels=0,
                    cmap='icefire',
                    # dendrogram_ratio=(0.2,0.001),
                    # row_cluster=False,col_cluster=False,
                    row_colors=[palette[5 + 3*int(i=='UC')] for i in tp],
                    figsize=(7, 70))
    plt.savefig('Quantification/ZExpressionHeatmap'+samp+'.png',dpi=936)