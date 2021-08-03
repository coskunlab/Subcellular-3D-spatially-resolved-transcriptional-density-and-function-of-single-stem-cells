#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:39:34 2021

@author: ryan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
import csv
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)


colors = ['skyblue','navy','yellow','cyan','magenta','green','lightpink','purple','gold','indigo','violet','crimson']
with open('RNAcountsBM.txt') as f:
    lines = f.readlines()
BM = [[int(j) for j in i.split()] for i in lines[1:]]
tpBM = ['BM' for i in BM]
with open('RNAcountsUC.txt') as f:
    lines = f.readlines()
UC = [[int(j) for j in i.split()] for i in lines[1:]]
tpUC = ['UC' for i in UC]
exp = BM+UC
tp=tpBM+tpUC
gene = lines[0].split()
gene = [i*(i!='ssp1')+'spp1'*(i=='ssp1') for i in gene]

for ctype,Dots in zip(['bm','uc'],[BM,UC]):
    dist = pd.DataFrame(
        np.log10(np.array(Dots)+0.1)[:,2:5],
        columns=gene[2:5])
    print(dist)
    dist.agg(['min', 'max', 'mean', 'std']).round(decimals=2)
    fig, ax = plt.subplots()
    dist.plot.kde(ax=ax, legend=False, color=colors[2:5])
    dist.plot.hist(density=True, bins=20, ax=ax, alpha=0.5, color=['gold','magenta','cyan'], title=ctype.upper()+' Expression: N = '+str(len(Dots)))
    ax.set_ylabel('Probability')
    ax.set_xlim((-0.99,4))
    ax.set_facecolor('white')
    fig=plt.gcf()
    fig.set_size_inches(7,5)
    plt.savefig('Quantification/'+ctype.upper()+'expressionDistributionLog.png')
    plt.show()
    
    d = pd.DataFrame([[k for i,k in zip(gene,j) if i in ('il8','il6','ccl11')] for j in Dots],columns=[i for i in gene if i in ('il8','il6','ccl11')])
    sns.violinplot(data=d/np.mean(d,axis=0))
    plt.title('C'+ctype)
    plt.xticks(rotation=45)
    fig=plt.gcf()
    fig.set_size_inches(5,5)
    plt.savefig('Quantification/'+ctype.upper()+'wineplot3D_cyto.png')
    plt.show()
    
palette = sns.color_palette()
for a,b,c in zip([BM,UC],[tpBM,tpUC],['BM','UC']):
    E = pd.DataFrame([i[2:5] for i in a],b,gene[2:5])
    sns.clustermap(E, metric="correlation",
                    xticklabels=1,yticklabels=0,
                    cmap='icefire',
                    row_cluster=True,col_cluster=False,
                    row_colors=[palette[5 + 3*int(i=='UC')] for i in b],
                    figsize=(7, 7))
    plt.savefig('Quantification/CytokineExpressionHeatmap'+c+'.png')
    


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