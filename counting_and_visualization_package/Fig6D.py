#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 03:41:06 2021

@author: ryan
"""


from os import system
_ = system('clear')
import os
from sklearn.cluster import KMeans
from os import path
from os import listdir
from os.path import isdir, isfile
import matplotlib.pyplot as plt
import numpy as np
import time
from PIL import Image, ImageDraw, ImageFont 
import seaborn as sns
from scipy import ndimage
import matplotlib
from skimage.morphology import local_maxima, disk, binary_dilation, binary_erosion
import itertools
import pandas as pd
from scipy.stats import zscore, pearsonr, linregress
from functools import reduce
from sklearn.manifold import TSNE

font_size_of_the_code = 12#24 #12
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)


Beginning = time.time()


dname = os.path.dirname(os.path.abspath(__file__)) #PYTHON SCRIPT DIRECTORY

"""
MAKE NEW FILE WITH PRETTY MUCH EVERYTHING DELETED

ALSO MAKE PREDICTION DATA ON NEW DATA
"""

os.chdir(dname)
e = 2.7182818284590452353602874713527
def translist(list_here):
    return [list(x) for x in zip(*list_here)]

def flatten(t):
    return [item for sublist in t for item in sublist]

def Measure(x):
    return x

os.chdir(dname)
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

for nor in [2]:#[0,1,2]:
    if nor==0: NORMIZE = lambda x,y: (np.log(x+0.1)/np.log(y+0.1)); nnnn = "log(x)/log(EEF2)"
    if nor==1: NORMIZE = lambda x,y: np.log( (x+0.1)/(y+0.1) ); nnnn = "log(x/EEF2)"
    if nor==2: NORMIZE = lambda x,y: (x+0.1)/(y+0.1); nnnn = "x/EEF2"
    # NORMIZE = lambda x,y: NORMIZE(x,y)**2
    makePicsOfARGMAX = 0
    makePicsOfBAYES = 1
    
    NormalizationGene='EEF2'#'GAPDH'#'ACTB'# EEF2 is the one selected from RNAseq
    Exp = {gene[i]: NORMIZE(np.array(exp[i]),np.array(Exp[NormalizationGene])) for i in range(len(gene))}
    def getProbFunc(geneSampled,ValSeen,EValSeen):
        LikelyhoodNotT = tuple(map(lambda x: abs(sum(map(lambda y: y - NORMIZE(ValSeen,EValSeen), x))), translist(Exp[geneSampled])))
        ProbT = tuple(map(lambda x: 0.01 + reduce(lambda x,y: x*y, [LikelyhoodNotT[i] for i in range(len(Time)) if (i != x)]), range(len(Time))))
        Porportionality = 1/sum(ProbT)
        return tuple(map(lambda x: x*Porportionality, ProbT))
    def getProbFuncBayes(geneSampled,ValSeen,EValSeen):
        GeneData = translist(Exp[geneSampled])
        mu = tuple(map(lambda x: sum(x)/len(x), GeneData))
        numtop = tuple(map(lambda x, m: sum(map(lambda y: (y - m)**2,x)), GeneData, mu))
        sig = tuple(map(lambda x, m: (x/(len(numtop)-1))**0.5 + 0.001, numtop, mu))
        ProbT = tuple(map(lambda m, s: 0.01 + np.nan_to_num((e**(-0.5*((m - NORMIZE(ValSeen,EValSeen))/s)**2)/s)), mu, sig))
        Porportionality = 1/sum(ProbT)
        return tuple(map(lambda x: x*Porportionality, ProbT))
        
    gcs2 = []
    cc2 = []
    gg2 = []
    for ff in ['006','034']:
        os.chdir('/home/ryan/Documents/back_sub/chondrocytes_count_data/count_data/'+ff+'/markers/')
        GENES = [f for f in listdir() if not ('_Mask' in f or 'nucleus' in f)]
        GENES.sort()
        msks = [j for j in listdir() if ('_Mask' in j)]
        n1 = 5
        n2 = 4
        n3 = 2
        n4 = 1
        convolvFilter = -disk(n1).astype(float)
        convolvFilter[n1-n2:len(convolvFilter)-n1+n2,n1-n2:len(convolvFilter)-n1+n2] += disk(n2).astype(float)
        convolvFilter[n1-n3:len(convolvFilter)-n1+n3,n1-n3:len(convolvFilter)-n1+n3] -= 0.1*disk(n3).astype(float)*np.sum(convolvFilter)/np.sum(disk(n3))
        convolvFilter[n1-n4:len(convolvFilter)-n1+n4,n1-n4:len(convolvFilter)-n1+n4] -= disk(n4).astype(float)*np.sum(convolvFilter)/np.sum(disk(n4))
        convolvFilter = convolvFilter/np.sum(abs(convolvFilter))
        cutoff = 1000
        
    
        def JustCount(file,mask):
            # I = tifffile.imread(file)
            # M = tifffile.imread(mask)
            I = Image.open(file)
            I = np.array(I)
            M = Image.open(mask)
            M = np.array(M)
            gene = file
            gene = gene[gene.find('_')+1:]
            gene = gene[:gene.find('_')]
            G = ndimage.convolve(I.astype(float),convolvFilter, mode='constant', cval=0)
            # print((np.max(G),np.mean(G),np.min(G)))
            B = G>1000
            # B = G/np.max(G) > 0.5
            DotCount = np.sum(local_maxima(B*G*M))
            return gene, DotCount#, G*(G>0)
        
        
        
        # gcs = [[JustCount(file,mask) for file in GENES[-2:-1]] for mask in msks[:1]]
        gcs = [[JustCount(file,mask) for file in GENES] for mask in msks]
        cc2.append([ff[1:]+'.'+i[4:i.find('_')] for i in msks])
        gcs2.append(gcs)
        gg2.append([i[i.find('_')+1:][:i[i.find('_')+1:].find('_')].upper() for i in GENES])
        gcs = [[i for i in j if i[0].upper() in gene] for j in gcs]
        # gcs = [[i for i in j if i[0].upper() in ('GAPDH','COL1A1','SOX9','RUNX1','EEF2','ACTB','SPP1',)] for j in gcs]
        # gcs = [[i for i in j if i[0].upper() in ('COL1A1','SOX9','EEF2','ACTB','SPP1',)] for j in gcs]
        # # I = gcs[0][0][-1]
        # # aaaaaaaaaaaaaaaaaa
        # # print([[i[:-1] for i in j if i[0].upper() == 'SOX9'] for j in gcs])
        # # print(gcs)
        # for i in gcs:
        #     plt.imshow(binary_dilation(i[0][-1]>1000,disk(20)))
        #     plt.show()
        # print([i[0] for i in gcs[0] if i[0].upper() in gene])#4
        # plt.bar(*translist(gcs[msks.index('cell4_Mask.tif')]))
        # plt.xticks(rotation = 45)
        # plt.show()
        
        os.chdir(dname)
        # print(translist(gcs)[8])
        exp = np.array([[i[1] for i in j] for j in gcs]).T
        
        
        
        
        
        
        sox9 = exp.T[translist(gcs[0])[0].index('sox9')]
        newexp = TSNE(n_components = 2).fit(exp).embedding_
        fig, ax = plt.subplots()
        plt.scatter(*newexp.T, c = sox9)
        # ax.autoscale(tight=True)
        ax.set_title('Color is sox9 expression')
        ax.axis("off")
        plt.savefig('Quantification/sox9TSNEchondrocytes2.png')
        plt.show()
        os.chdir('/home/ryan/Documents/back_sub/chondrocytes_count_data/count_data/'+ff+'/markers/')
        
        for makePicsOfBAYES in (0,1):
            for makePicsOfARGMAX in (0,1):
                def JustPredict(gc):
                    geens,cts = [i[0].upper() for i in gc],[i[1] for i in gc]
                    if makePicsOfBAYES:
                        getProbFuncR = getProbFuncBayes
                    else:
                        getProbFuncR = getProbFunc
                    Pthis = [getProbFuncR(ge,ctss,cts[geens.index(NormalizationGene)]) for ge, ctss in zip(geens,cts) if ge != 'EEF2']
                    # print(Pthis)
                    Ptot = np.prod(np.array(Pthis),axis=0)
                    if makePicsOfARGMAX:
                        t = Time[np.argmax(Ptot)]
                    else:
                        t = round(np.ceil(np.sum(Time*Ptot/np.sum(Ptot))))
                    return t
                times = [JustPredict(gc) for gc in gcs]
                # print(translist([times] + translist(translist(gcs)[8])))
                print((times, 'bayes', makePicsOfBAYES, 'argmax', makePicsOfARGMAX, 'model',nnnn))
            plt.bar([0,1],[np.mean([i for i,j in zip(times,sox9) if j<100]),np.mean([i for i,j in zip(times,sox9) if j>=100])])
            plt.title('bayes: '+str(makePicsOfBAYES)+', argmax: '+str(makePicsOfARGMAX)+', '+nnnn)
            os.chdir(dname)
            plt.savefig('Quantification/SOX9cell'+ff+'normtype'+str(len(nnnn))+'bayes'+str(makePicsOfBAYES)+'.png')
            plt.show()
            os.chdir('/home/ryan/Documents/back_sub/chondrocytes_count_data/count_data/'+ff+'/markers/')
        
        # translist(gcs)[-2]
        
        
        
        # plt.bar([0,1],[np.mean([i for i,j in zip(times,sox9) if j<100]),np.mean([i for i,j in zip(times,sox9) if j>=100])])
        # plt.show()
    os.chdir(dname)
    
    
    
    exp = zscore(np.array([[i[1] for i in j] for j in flatten(gcs2)]).T,axis=1)
    gg = gg2[0]
    E = pd.DataFrame(exp,gg,flatten(cc2))
    g = sns.clustermap(E, metric="correlation",
                    xticklabels=1,yticklabels=1,
                    cmap='icefire',
                    dendrogram_ratio=(0.2,0.001),
                    cbar_pos=(0,0.1,0.05,0.9),
                    cbar_kws={'label': ""},
                    # row_cluster=False,col_cluster=False,
                    # figsize=(8, 10),
                    )
    g.ax_heatmap.set_title('RNA expression',pad=0.1)
    g.figsize=(10,10)
    plt.savefig('Quantification/newChondrocytes.png')
    
    




Runtime = time.time() - Beginning
print('Runtime: '+str(Runtime//60)[:-2]+':'+'0'*(2-len(str(Runtime%60 //1)[:-2]))+str(Runtime%60 //1)[:-2])