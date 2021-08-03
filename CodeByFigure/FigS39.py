#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:16:06 2021

@author: ryan
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
# import umap.umap_ as umap
import umap
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)

def translist(list_here):
    return [list(x) for x in zip(*list_here)]


def flatten(t):
    return [item for sublist in t for item in sublist]

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
    totdata = pd.DataFrame(AllStuff,index=['BM','UC','fc','BM std','UC std'],columns=gene).T
    """
    UMAP ANALYSIS OF THE HCR CELL DATA
    """
    AllData = np.array(flatten([BM,UC]))
    reducer = umap.UMAP(random_state=2, n_neighbors=20, metric='correlation')
    AllData = StandardScaler().fit_transform(AllData)
    embedding = reducer.fit_transform(AllData)
    plt.scatter(*translist(embedding),c = [0]*len(BM)+[1]*len(UC), cmap='RdBu')
    plt.title(samp)
    plt.savefig('Quantification/UMAPcells'+samp+'.png')
    plt.show()
    
    if samp == '':
        e1 = np.mean(np.array(flatten([BM,UC]))[np.array(embedding).T[0]>np.mean(np.array(embedding).T[0])],axis=0)
        e2 = np.mean(np.array(flatten([BM,UC]))[np.array(embedding).T[0]<=np.mean(np.array(embedding).T[0])],axis=0)
        print(np.array(gene)[(e1-e2)>0])
        print(np.array(gene)[(e1-e2)<=0])
        genechanges = tuple(zip(gene,[(str(i)[:5]) for i in 100*(e1-e2)/np.mean(np.array(flatten([BM,UC])),axis=0)]))
        print('\n'.join([': '.join(i) for i in genechanges]))
        print('\n'.join([k+': '+str(int(i*10)/10)+'    |    '+str(int(j*10)/10) for k,i,j in zip(gene,e1,e2)]))
        
    # inds = KMeans(n_clusters=2, random_state=0).fit(embedding).labels_
l2fc = list(map(lambda a: [np.log2(y/x) for x,i in zip(*a[0]) for y,j in zip(*a[1]) if i == j], mm))
gene = [i for x,i in zip(*mm[0][0]) for y,j in zip(*mm[0][1]) if i == j]
l2fc = [[str(i)[:str(float(i)).find('.')+3] for i in j] for j in l2fc]
totdata = pd.DataFrame((l2fc),index=['fcBM','fcUC'],columns=gene).T
