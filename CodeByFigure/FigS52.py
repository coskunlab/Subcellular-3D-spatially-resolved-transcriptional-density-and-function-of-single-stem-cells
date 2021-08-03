#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:44:58 2021

@author: ryan
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)

colors = ['skyblue','navy','yellow','cyan','magenta','green','lightpink','purple','gold','indigo','violet','crimson']

for dDdate,samp in zip(['05_12','05_12'],['low_dense','high_dense']):
    with open('RNAcounts'+dDdate+'BM'+samp+'.txt') as f:
        lines = f.readlines()
    BM = [[float(j) for j in i.split()] for i in lines[1:]]
    tpBM = ['BM' for i in BM]
    with open('RNAcounts'+dDdate+'UC'+samp+'.txt') as f:
        lines = f.readlines()
    UC = [[float(j) for j in i.split()] for i in lines[1:]]
    tpUC = ['UC' for i in UC]
    exp = BM+UC
    tp=tpBM+tpUC
    gene = lines[0].split()
    
    
    for ctype,Dots in zip(['bm','uc'],[BM,UC]):
        dist = pd.DataFrame(
            np.log10(np.array(Dots)+0.1)[:,3:6],
            columns=gene[3:6])
        dist.agg(['min', 'max', 'mean', 'std']).round(decimals=2)
        fig, ax = plt.subplots()
        dist.plot.kde(ax=ax, legend=False, color=['gold','magenta','cyan'])
        dist.plot.hist(density=True, bins=20, ax=ax, alpha=0.5, color=colors[2:5], title=samp+ctype.upper()+' Expression: N = '+str(len(Dots)))
        ax.set_ylabel('Probability')
        ax.set_xlim((-0.99,4))
        ax.set_facecolor('white')
        fig=plt.gcf()
        fig.set_size_inches(7,5)
        plt.savefig('Quantification/'+ctype.upper()+'expressionDistribution'+dDdate+samp+'.png')
        plt.show()


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
    gene = [i*(i!='ssp1')+'spp1'*(i=='ssp1') for i in gene]
    
    for ctype,Dots in zip(['bm','uc'],[BM,UC]):
        if samp == '':
            dist = pd.DataFrame(
                np.log10(np.array(Dots)+0.1)[:,2:5],
                columns=gene[2:5])
        else:
            dist = pd.DataFrame(
                np.hstack((np.log10(np.array(Dots)+0.1)[:,3:5],np.log10(np.array(Dots)+0.1)[:,:1])),
                columns=gene[3:5]+gene[:1])
        print(dist)
        dist.agg(['min', 'max', 'mean', 'std']).round(decimals=2)
        fig, ax = plt.subplots()
        dist.plot.kde(ax=ax, legend=False, color=colors[2:5])
        dist.plot.hist(density=True, bins=20, ax=ax, alpha=0.5, color=['gold','magenta','cyan'], title=samp+ctype.upper()+' Expression: N = '+str(len(Dots)))
        ax.set_ylabel('Probability')
        ax.set_xlim((-0.99,4))
        ax.set_facecolor('white')
        fig=plt.gcf()
        fig.set_size_inches(7,5)
        plt.savefig('Quantification/'+ctype.upper()+'expressionDistributionLog'+samp+'.png')
        plt.show()