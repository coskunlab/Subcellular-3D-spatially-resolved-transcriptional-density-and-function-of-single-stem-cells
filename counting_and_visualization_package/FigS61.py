#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:55:50 2021

@author: ryan
"""
import numpy as np
import pandas as pd
import umap
# import umap.umap_ as umap
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)





def flatten(x):
    return [i for j in x for i in j]

def translist(x):
    return [list(i) for i in zip(*x)]

"""
CREATES THE SUPPLIMENTARY TABLES WITH FOLD CHANGE UNDER totdata

AND UMAP ANALYSIS OF EXPRESSSION OVER CELLS
"""
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



"""
UMAP of Chondrogenesis from RNA-Seq
"""
genesRaw = [i for i in open('PredictionData/ChondrogenesisOnlineData.txt',"r")]

donor = [0,1,2]
Time = np.array([0,1,3,7,14,21])
TIME = [i for i,j in enumerate(Time)]*3
cats = genesRaw[0].split()
geneind = cats.index('geneName')
gene = [i.split()[geneind] for i in genesRaw[1:]]
timeind = translist([[cats.index(i)-2 for i in cats if 'D'+str(j)+'D' in i] for j in [1,2,3]])

exp = [i.split()[2:] for i in genesRaw[1:]]
exp = tuple(map(lambda  x: [[float(x[k]) for k in l] for l in translist(timeind)], exp))
Exp = {gene[i]: exp[i] for i in range(len(gene))} 

GENESFROMDATA= [i for i in """ACTB
    GAPDH
    EEF2
    NANOG
    SOX9
    RUNX1
    SPP1
    COL1A1
    COL5A2
    IL6
    IL8
    CCL11
    CXCR4
    CD274
    MKI67
    MALAT1
    """.split() if i in gene]
HCRexp = np.array([flatten(i) for i in [Exp[i] for i in GENESFROMDATA]]).T

np.random.seed(seed=0)
reducer = umap.UMAP(random_state=np.random.randint(2**8), n_components=2, n_neighbors=len(TIME)-1, min_dist=0.5, metric='correlation')
AllData = StandardScaler().fit_transform(HCRexp)
embedding = reducer.fit_transform(AllData)
fig, ax = plt.subplots()
ax.scatter(*translist(embedding),c = TIME, cmap='gnuplot')
fig.suptitle('UMAP of RNA-Seq: '+' '.join(GENESFROMDATA))

fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=min(TIME), vmax=max(TIME)), cmap='gnuplot'), ax=ax)
plt.show()

NGenes = 12
NTests = (3,3,)
fig, ax = plt.subplots(*NTests)
for n in range(np.prod(NTests)):
    EXP = np.random.permutation([flatten(i) for i in exp])[:NGenes].T
    AllData = StandardScaler().fit_transform(EXP)
    embedding = reducer.fit_transform(AllData)
    ax[n%NTests[0],int(np.floor(n/NTests[0]))].scatter(*translist(embedding),c = TIME, cmap='gnuplot')#, cmap='coolwarm')#, cmap='icefire')
    ax[n%NTests[0],int(np.floor(n/NTests[0]))].set_xticks([])
    ax[n%NTests[0],int(np.floor(n/NTests[0]))].set_yticks([])
fig.suptitle('UMAP of RNA-Seq: '+str(NGenes)+' genes')
plt.show()# Bright is Stem, Black is Chondrocyte

