#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:45:21 2021

@author: ryan
"""
import numpy as np
import pandas as pd
import csv





"""
TABLE 2,3,4,5

2,3,4 are saves in totdata

5 is saved in tble
"""
def flatten(x):
    return [i for j in x for i in j]

totdata = []
ss, mm = [[[[]],[[]]],[[[]],[[]]]], [[[[]],[[]]],[[[]],[[]]]]
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
    s, m = [[],[]], [[],[]]
    for n, Dots in enumerate([BM,UC]):
        # ss[n][(samp=='3D')] = np.std([[k for k in (j)] for j in Dots],axis=0)
        mm[n][(samp=='3D')] = (np.mean([[k for k in (j)] for j in Dots],axis=0), gene)
        s[n] = np.std([[k for k in (j)] for j in Dots],axis=0)
        m[n] = np.mean([[k for k in (j)] for j in Dots],axis=0)
    l2fc = list(map(lambda x,y: np.log2(y/x), *m))
    AllStuff = m+[l2fc]+s
    AllStuff = [[str(i)[:str(float(i)).find('.')+3] for i in j] for j in AllStuff]
    totdata.append(pd.DataFrame(AllStuff,index=['BM','UC','fc','BM std','UC std'],columns=gene).T)
l2fc = list(map(lambda a: [np.log2(y/x) for x,i in zip(*a[0]) for y,j in zip(*a[1]) if i == j], mm))
gene = [i for x,i in zip(*mm[0][0]) for y,j in zip(*mm[0][1]) if i == j]
l2fc = [[str(i)[:str(float(i)).find('.')+3] for i in j] for j in l2fc]
totdata.append(pd.DataFrame((l2fc),index=['fcBM','fcUC'],columns=gene).T)


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

UCind = [i for i,j in enumerate(Sample) if 'UC' in j.upper()]
BMind = [i for i,j in enumerate(Sample) if 'BM' in j.upper()]
UCamnt = np.mean(np.array(Amount)[min(UCind):max(UCind)+1,:],axis=0)
BMamnt = np.mean(np.array(Amount)[min(BMind):max(BMind)+1,:],axis=0)
tble = [[round(i) for i in UCamnt], [round(i) for i in BMamnt],[round(i) for i in Zeros]]
tble = pd.DataFrame(tble, columns = Cytokine, index = ['UC','BM','Negative Control']).T
