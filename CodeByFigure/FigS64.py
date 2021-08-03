#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:19:48 2021

@author: ryan
"""


import os
from os import path
from os import listdir
from os.path import isdir
import matplotlib.pyplot as plt
import numpy as np
import time
from PIL import Image, ImageDraw, ImageFont 
import seaborn as sns
from scipy import ndimage
import matplotlib
from functools import reduce

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
matplotlib.rc('font', **font)


Beginning = time.time()


dname = os.path.dirname(os.path.abspath(__file__)) #PYTHON SCRIPT DIRECTORY



os.chdir(dname)
e = 2.7182818284590452353602874713527
def translist(list_here):
    return [list(x) for x in zip(*list_here)]

def flatten(t):
    return [item for sublist in t for item in sublist]

def Measure(x):
    return x

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109503
genesRaw = [i for i in open('PredictionData/ChondrogenesisOnlineData.txt',"r")]

donor = [0,1,2]
Time = np.array([0,1,3,7,14,21])

cats = genesRaw[0].split()
geneind = cats.index('geneName')
gene = [i.split()[geneind] for i in genesRaw[1:]]
timeind = translist([[cats.index(i)-2 for i in cats if 'D'+str(j)+'D' in i] for j in [1,2,3]])

exp = [i.split()[2:] for i in genesRaw[1:]]
exp = tuple(map(lambda  x: [[Measure(float(x[k])) for k in l] for l in translist(timeind)], exp))

Exp = {gene[i]: exp[i] for i in range(len(gene))}




GenesInAllHCR = ['GAPDH', 'ACTB', 'IL8', 'IL6', 'CCL11', 'COL1A1', 'NANOG', 'SOX9', 'EEF2', 'SSP1', 'RUNX1', 'PDL1']
GenesInAllHCR = [i for i in GenesInAllHCR if i in ['PDL1', 'SOX9', 'CL1A1', 'ACTB', 'IL8', 'IL6', 'CCL11', 'RUNX1', 'SPP1', 'EEF2']]

GenesInAllHCR = [i for i in GenesInAllHCR if i in gene]
# print(sorted([(np.mean(Exp[i]),i) for i in GenesInAllHCR],reverse=True))
# print(sorted([(np.std(Exp[i])/np.mean(Exp[i]),i) for i in GenesInAllHCR]))
# NormGene = sorted([(np.mean(Exp[i]),i) for i in GenesInAllHCR],reverse=True)[0][-1]
NormGene = sorted([(np.std(Exp[i])/np.mean(Exp[i]),i) for i in GenesInAllHCR])[0][-1]
NormNum = np.array(Exp[NormGene])

for which in (-1,0,1,):
    if which == -1:
        normize = lambda x,y: (x+0.1)/(y+0.1)
        nm=''
    elif which == 0:
        normize = lambda x,y: np.log(x+0.1)/np.log(y+0.1)
        nm='dl'
    elif which == 1:
        normize = lambda x,y: np.log( (x+0.1)/(y+0.1) )
        nm='ld'
    
    Exp = {gene[i]: normize(np.array(exp[i]),NormNum) for i in range(len(gene)) if not gene[i] == NormGene}
    gene = [i for i in gene if not i==NormGene]
    
    
    def getProbFunc(geneSampled,ValSeen):
        LikelyhoodNotT = tuple(map(lambda x: abs(sum(map(lambda y: y - ValSeen, x))), translist(Exp[geneSampled])))
        ProbT = tuple(map(lambda x: 0.01 + reduce(lambda x,y: x*y, [LikelyhoodNotT[i] for i in range(len(Time)) if (i != x)]), range(len(Time))))
        Porportionality = 1/sum(ProbT)
        return tuple(map(lambda x: x*Porportionality, ProbT))
    VarianceFromNoise = 0
    def getProbFuncBayes(geneSampled,ValSeen):
        GeneData = translist(Exp[geneSampled])
        mu = tuple(map(lambda x: sum(x)/len(x), GeneData))
        numtop = tuple(map(lambda x, m: sum(map(lambda y: (y - m)**2,x)), GeneData, mu))
        sig = tuple(map(lambda x, m: (x/(len(numtop)-1) + VarianceFromNoise)**0.5, numtop, mu))
        ProbT = tuple(map(lambda m, s: 0.01 + e**(-0.5*((m - ValSeen)/s)**2)/s, mu, sig))
        # ProbT = tuple(map(lambda m, s: e**(-0.5*((m - ValSeen)/s)**2)/s, mu, sig))
        Porportionality = 1/sum(ProbT)
        return tuple(map(lambda x: x*Porportionality, ProbT))

    
    FinaleSample = ['col1a1'.upper(), 106,
        'eef2'.upper(), 804,
        'spp1'.upper(), 68,
        'il8'.upper(), 62,
        'il6'.upper(), 48,
        'ccl11'.upper(), 37,
        'runx1'.upper(), 156,
        'pdl1'.upper(), 3,
        'sox9'.upper(), 62,#248,#
        'actb'.upper(), 45,#89,#
        ]
    FinaleSample = tuple(zip(FinaleSample[::2],FinaleSample[1::2]))
    FinaleSample = [(i[0],normize(i[1],FinaleSample[1][1])) for i in FinaleSample if (not (i==FinaleSample[1])) and (i[0] in gene)]
    
    print('\n\n\n')
    for pfunc in [getProbFunc,getProbFuncBayes]:
        FinalExpectation = np.sum(Time*np.prod([pfunc(g,d) for g,d in FinaleSample],axis=0)/np.sum(np.prod([pfunc(g,d) for g,d in FinaleSample],axis=0)))
        FinalProbs = np.prod([pfunc(g,d) for g,d in FinaleSample],axis=0)/np.sum(np.prod([pfunc(g,d) for g,d in FinaleSample],axis=0))
        print(FinalExpectation)
    
    print('Normalization based on: '+NormGene)
    print('\n')
    palette = sns.color_palette()
    for cult in ['','3D']:
        with open('RNAcountsBM'+cult+'.txt') as f:
            lines = f.readlines()
        BM = [[float(j) for j in i.split()] for i in lines[1:]]
        tpBM = ['BM' for i in BM]
        with open('RNAcountsUC'+cult+'.txt') as f:
            lines = f.readlines()
        UC = [[float(j) for j in i.split()] for i in lines[1:]]
        tpUC = ['UC' for i in UC]
        
        HCRgene = [i.upper() for i in lines[0].split()]
        NormGeneInd = HCRgene.index(NormGene)
        UC = [[normize(j,i[NormGeneInd]) for j in i] for i in UC if not i[NormGeneInd]==0]
        BM = [[normize(j,i[NormGeneInd]) for j in i] for i in BM if not i[NormGeneInd]==0]
        
        BM = [ translist(tuple(filter(lambda x: x[0] in gene, zip(HCRgene,i))))[1] for i in BM]
        UC = [ translist(tuple(filter(lambda x: x[0] in gene, zip(HCRgene,i))))[1] for i in UC]
        HCRgene = tuple(filter(lambda x: x in gene, HCRgene))
        
        
        
        
        
        BMtime = [ np.sum(Time*np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in BM]
        UCtime = [ np.sum(Time*np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in UC]
        bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
        uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
        x=x[:-1]
        ax = plt.subplot(111)
        line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
        line2=ax.bar(x, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
        ax.autoscale(tight=True)
        ax.set_title('Prediction MLP '+cult)
        ax.set_xlabel('Differentiation Time')
        ax.legend([line1, line2], ['BM', 'UC'])
        plt.savefig('Quantification/NormDifferentiaionPrediction'+nm+'MLP'+cult+'.png')
        plt.show()
        
        
        BMtime = [ Time[np.argmax(np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0))] for c in BM]
        UCtime = [ Time[np.argmax(np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0))] for c in UC]
        bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
        uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
        x=x[:-1]
        ax = plt.subplot(111)
        line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
        line2=ax.bar(x, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
        ax.autoscale(tight=True)
        ax.set_title('Prediction MLP '+cult)
        ax.set_xlabel('Differentiation Time')
        ax.legend([line1, line2], ['BM', 'UC'])
        plt.savefig('Quantification/NormDifferentiaionPrediction'+nm+'MLP'+cult+'ARGMAX.png')
        plt.show()
        
        BMtime = [ np.sum(Time*np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in BM]
        UCtime = [ np.sum(Time*np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in UC]
        bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
        uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
        x=x[:-1]
        ax = plt.subplot(111)
        line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
        line2=ax.bar(x+0.2, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
        ax.autoscale(tight=True)
        ax.set_title('Prediction Bayes '+cult)
        ax.set_xlabel('Differentiation Time')
        ax.legend([line1, line2], ['BM', 'UC'])
        plt.savefig('Quantification/NormDifferentiaionPrediction'+nm+'Bayes'+cult+'.png')
        plt.show()
        
        BMtime = [ Time[np.argmax(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))] for c in BM]
        UCtime = [ Time[np.argmax(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))] for c in UC]
        bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
        uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
        x=x[:-1]
        ax = plt.subplot(111)
        line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
        line2=ax.bar(x+0.2, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
        ax.autoscale(tight=True)
        ax.set_title('Prediction Bayes '+cult)
        ax.set_xlabel('Differentiation Time')
        ax.legend([line1, line2], ['BM', 'UC'])
        plt.savefig('Quantification/NormDifferentiaionPrediction'+nm+'Bayes'+cult+'ARGMAX.png')
        plt.show()
        
        print(cult+' Predictions From: '+' '.join(HCRgene))





Runtime = time.time() - Beginning
print('Runtime: '+str(Runtime//60)[:-2]+':'+'0'*(2-len(str(Runtime%60 //1)[:-2]))+str(Runtime%60 //1)[:-2])