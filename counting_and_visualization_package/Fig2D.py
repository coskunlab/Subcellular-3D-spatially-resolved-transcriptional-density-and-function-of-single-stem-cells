#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:04:02 2021

@author: ryan
"""
import numpy as np
import matplotlib.pyplot as plt
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

for ctype,Dots in zip(['bm','uc'],[BM,UC]):
    plt.bar(gene,np.mean(Dots,axis=0), yerr=np.std(Dots,axis=0), color=colors)
    plt.title(ctype.upper()+' mean: N = '+str(len(Dots)))
    plt.xticks(rotation=45)
    fig=plt.gcf()
    fig.set_size_inches(6,5)
    plt.savefig('Quantification/'+ctype.upper()+'mean.png')
    plt.show()