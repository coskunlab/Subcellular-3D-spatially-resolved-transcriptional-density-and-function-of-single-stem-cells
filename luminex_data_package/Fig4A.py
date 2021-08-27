#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:28:49 2021

@author: ryan
"""
import numpy as np
from scipy.stats import zscore
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)



csv_file = open('PredictionData/AllMedianData.csv')
read_tsv = csv.reader(csv_file, delimiter=",")
data = []
for row in read_tsv:
    data.append(row)
csv_file.close()

Cytokine = data[0][1:]
Sample = tuple(map(lambda x: x[0], data))[1:]
Sample = Sample[:-2]
Amount = tuple(map(lambda x: tuple(map(lambda y: float(y),x[1:])), data[1:]))
Zeros = tuple(map(lambda x, y: (x+y)*0.5, Amount[-2], Amount[-1]))

Amount = tuple(map(lambda x: tuple(map(lambda y, z: max([y-z,0]), x, Zeros)), Amount[:-2]))
Sample = [Sample[i] for i in range(len(Sample)) if not i==10]
Amount = [[j - k for j,k in zip(Amount[i],Amount[10])] for i in range(len(Amount)) if not i==10]

def ModifyName(x):
    p1 = 'UC'*('UC' in x.upper())+'BM'*('BM' in x.upper())
    p2 = x[x.find('P'):x.find('P')+2]
    p3 = 'LD'*('LOW DENSITY' in x.upper())+'HD'*('HIGH DENSITY' in x.upper())
    if not p3:
        p3 = 'VD'
    return p1+'-'+p2+'-'+p3
samps = [ModifyName(i) for i in Sample]

for m,l,t in zip([lambda x: np.log10(np.array(x)*(np.array(x)>0)+0.1), lambda x:np.array(zscore(x,axis=0))],['log10(AFU)','Amount'],['','Zscored']):
    E = pd.DataFrame(m(Amount).T,Cytokine,samps)
    g = sns.clustermap(E, metric="correlation",
                    xticklabels=1,yticklabels=1,
                    cmap='icefire',
                    dendrogram_ratio=(0.2,0.001),
                    cbar_pos=(0,0.1,0.05,0.9),
                    cbar_kws={'label': l},
                    # row_cluster=False,col_cluster=False,
                    # figsize=(8, 10),
                    )
    g.ax_heatmap.set_title('Luminex Data '+t,pad=0.1)
    g.figsize=(10,10)
    plt.savefig('Quantification/LuminexTotal'+t+'.png')
