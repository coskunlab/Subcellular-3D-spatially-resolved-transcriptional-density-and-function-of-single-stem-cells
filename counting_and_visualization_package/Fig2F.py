#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:18:42 2021

@author: ryan
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)



def flatten(x):
    return [i for j in x for i in j]

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



for gn in gene:
    fig, ax = plt.subplots()
    d = pd.DataFrame([i+['BM'] for i in BM]+[i+['UC'] for i in UC],columns=gene+['Type'])
    d1 = [[i,'BM',gene[g%len(gene)]] for g,i in enumerate(flatten(BM)) if gene[g%len(gene)] == gn]
    d2 = [[i,'UC',gene[g%len(gene)]] for g,i in enumerate(flatten(UC)) if gene[g%len(gene)] == gn]
    d = pd.DataFrame(d1+d2,columns=['Amount','Type','Gene'])
    sns.violinplot(x=d.Amount,y=d.Gene,hue=d.Type,split=True)
    plt.xticks(rotation=45)
    fig=plt.gcf()
    fig.set_size_inches(7,5)
    plt.savefig('Quantification/Wineplot'+gn.upper()+'.png')
    plt.show()