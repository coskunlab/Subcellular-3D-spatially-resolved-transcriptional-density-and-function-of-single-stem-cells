#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:20:51 2021

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



def flatten(t):
    return [item for sublist in t for item in sublist]

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


fig, ax = plt.subplots()
d = pd.DataFrame([i+['BM'] for i in BM]+[i+['UC'] for i in UC],columns=gene+['Type'])
d1 = [[i/np.mean(d,axis = 0)[gene[g%len(gene)]],'BM',gene[g%len(gene)]] for g,i in enumerate(flatten(BM))]
d2 = [[i/np.mean(d,axis = 0)[gene[g%len(gene)]],'UC',gene[g%len(gene)]] for g,i in enumerate(flatten(UC))]
d = pd.DataFrame(d1+d2,columns=['Amount','Type','Gene'])
sns.violinplot(x=d.Amount,y=d.Gene,hue=d.Type,split=True)
plt.xlim((0,15))
plt.xticks(ticks=[])
fig=plt.gcf()
fig.set_size_inches(6.5,5)
plt.savefig('Quantification/Wineplot.png')
plt.show()