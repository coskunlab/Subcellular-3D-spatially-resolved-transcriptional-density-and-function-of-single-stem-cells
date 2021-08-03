#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:48:57 2021

@author: ryan
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import zscore, pearsonr, linregress
import csv
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)


def translist(x):
    return [list(i) for i in zip(*x)]




    

"""
CYTOKINE DATA
"""
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

  
def UCsampleHCRfinder(x):
      return max([x[0].endswith(', alpha-mem, collagen, '+i) and not 'P2' in x[0] and not 'BM' in x[0] for i in ['low density 1','low density 2','high density']])
def BMsampleHCRfinder(x):
      return max([x[0].endswith(', alpha-mem, collagen, '+i) and not 'P2' in x[0] and not 'UC' in x[0] for i in ['low density 1','low density 2','high density 1'] if not 'P2' in i])
Samps, amount = tuple(translist(sorted(translist([Sample, Amount]), key = lambda x: 'uc' in x[0].lower())))
UCSamps, UCamount = tuple(translist(filter(UCsampleHCRfinder, translist([Samps, amount]))))
BMSamps, BMamount = tuple(translist(filter(BMsampleHCRfinder, translist([Samps, amount]))))
Samps = UCSamps + BMSamps
amount = UCamount + BMamount
HCRcytes, amount = tuple(translist(filter(lambda x: x[0] in ['Eotaxin','IL-8','IL-6'],translist([Cytokine, translist(amount)]))))
def FixNames(x):
    return x[:2]+' '+x[x.rfind(',')+2].upper()+'D'+x[-1]*(x[x.rfind(',')+2]=='l')
Samps = [FixNames(i) for i in Samps]
HCRcytes[0] = HCRcytes[0]+'\n/ccl11'
amount = [i[::-1] for i in amount]
Samps = Samps[::-1]
HCRcytes = HCRcytes[::-1]
amount = amount[::-1]
y = amount


totsexp = []
for typ in ['BM','UC']:
    for samp in ['05_12'+typ+'low_dense',typ,'05_12'+typ+'high_dense']:
        with open('RNAcounts'+samp+'.txt') as f:
            lines = f.readlines()
        gene = lines[0].split()
        eee = np.mean([[float(j) for j in i.split()] for i in lines[1:]],axis=0)
        totsexp.append([eee[gene.index('il8')],eee[gene.index('il6')],eee[gene.index('ccl11')]])    
amount = np.array(totsexp).T
x = amount

"""
Calculates a Pearson correlation coefficient and the p-value for testing non-correlation.

The Pearson correlation coefficient measures the linear relationship between two datasets.
Strictly speaking, Pearson’s correlation requires that each dataset be normally distributed.
Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation.
Correlations of -1 or +1 imply an exact linear relationship.
Positive correlations imply that as x increases, so does y.
Negative correlations imply that as x increases, y decreases.

The p-value roughly indicates the probability of an uncorrelated system producing datasets that have a Pearson correlation at least as extreme as the one computed from these datasets.
The p-values are not entirely reliable but are probably reasonable for datasets larger than 500 or so.

THIS FIGURE SHOWS HOW WELL HCR CORROLATES WITH LUMINEX
"""
fig, ax = plt.subplots(figsize=(8, 6))
cCoOlLoOrRsS = ['yellow','cyan','magenta']
leg = []
for i,j in enumerate(cCoOlLoOrRsS):
    plt.scatter(x[i]/np.max(x[i]),y[i]/np.max(y[i]),c=j,linewidths=20)
    leg.append(HCRcytes[i]+': r = '+str(pearsonr(x[i],y[i])[0])[:5]+', p = '+str(pearsonr(x[i],y[i])[1])[:5])
plt.legend(leg)
plt.xlabel('HCR')
plt.ylabel('Luminex')
plt.title('Pearson Correlation')
ax.set_xticks([])
ax.set_yticks([])
plt.savefig('Quantification/HCR_Luminex_Pearson.png')
plt.show()

