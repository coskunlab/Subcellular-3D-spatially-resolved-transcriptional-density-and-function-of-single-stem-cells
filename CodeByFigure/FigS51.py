#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:42:44 2021

@author: ryan
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)





colors = ['skyblue','navy','yellow','cyan','magenta','green','lightpink','purple','gold','indigo','violet','crimson']
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
                np.array(Dots)[:,2:5],
                columns=gene[2:5])
        else:
            dist = pd.DataFrame(
                np.hstack((np.array(Dots)[:,3:5],np.array(Dots)[:,:1])),
                columns=gene[3:5]+gene[:1])
        print(dist)

        plt.bar([i for i in gene if i in ('il8','il6','ccl11')],np.mean([[k for i,k in zip(gene,j) if i in ('il8','il6','ccl11')] for j in Dots],axis=0), yerr=np.std([[k for i,k in zip(gene,j) if i in ('il8','il6','ccl11')] for j in Dots],axis=0)/len(Dots), color=[i for i,j in zip(colors,gene) if j in ('il8','il6','ccl11')])
        plt.title(samp+ctype.upper()+' mean: N = '+str(len(Dots)))
        plt.xticks(rotation=45)
        fig=plt.gcf()
        fig.set_size_inches(6,5)
        plt.savefig('Quantification/'+ctype.upper()+'mean'+samp+'_cyto.png')
        plt.show()
        
        d = pd.DataFrame([[k for i,k in zip(gene,j) if i in ('il8','il6','ccl11')] for j in Dots],columns=[i for i in gene if i in ('il8','il6','ccl11')])
        sns.violinplot(data=d/np.mean(d,axis=0))
        plt.title(samp+ctype)
        plt.xticks(rotation=45)
        fig=plt.gcf()
        fig.set_size_inches(5,5)
        plt.savefig('Quantification/'+ctype.upper()+'wineplot'+samp+'_cyto.png')
        plt.show()
