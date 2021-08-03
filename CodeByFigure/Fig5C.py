#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:06:46 2021

@author: ryan
"""

import os
from os import path
from os import listdir
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont 
from scipy import ndimage
from skimage.morphology import local_maxima, disk
from functools import reduce
import colorsys
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)



dname = os.path.dirname(os.path.abspath(__file__)) #PYTHON SCRIPT DIRECTORY


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


"""
HOW THE GENE STUFF IS NORMALIZED FOR ENCLUDING BOTH DATASETS
"""
# NORMIZE = lambda x,y: (np.log(x+0.1)/np.log(y+0.1))
NORMIZE = lambda x,y: np.log( (x+0.1)/(y+0.1) )
# NORMIZE = lambda x,y: (x+0.1)/(y+0.1)
makePicsOfARGMAX = 0
makePicsOfBAYES = 0



NormalizationGene='EEF2'
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
    ProbT = tuple(map(lambda m, s: 0.01 + np.nan_to_num((e**(-0.5*((m - NORMIZE(ValSeen,EValSeen))/s)**2)/s),nan=0), mu, sig))
    Porportionality = 1/sum(ProbT)
    return tuple(map(lambda x: x*Porportionality, ProbT))


FOLDER = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/12_gene_data/'
FOLDER3D = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/3D dataset/'
FOLDER7 = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/7_gene_data/'
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


def makePretty(I):
    I = np.array(I,dtype=float)
    I = (I-np.min(I))/(np.max(I)-np.min(I))*255
    return np.reshape(I, (I.shape[0],I.shape[1],1))


def makeImage7(mask, ims, genes, fname):
    M = plt.imread(mask)
    rgb = 0
    brightness = 1
    RNAcount = []
    for cc, file in enumerate(ims):
        I = plt.imread(file).astype('float64')
        
        G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
        dd = makePretty(G*(G>0)/np.max(G)*M)
        HSV = colorsys.hsv_to_rgb(cc/len(ims), 1, brightness)
        rgb = rgb + np.concatenate( [dd*HSV[0],dd*HSV[1],dd*HSV[2]],axis=2)
        B = G>cutoff
        RNAcount.append(np.sum(local_maxima(B*G*M)))
    img = Image.fromarray(np.clip(rgb*5,0,255).astype(np.uint8))
    imWidth, imHeight = img.size
    draw = ImageDraw.Draw(img) #2304 pixel = 249.6 micron
    l = (imWidth / 20, imHeight - imHeight / 20,
          imWidth / 20 + round(50*2304/249.6), imHeight - imHeight / 20)
    draw.line(l, fill='white', width=10)
    text = '50 \u03BCm'
    fontsize = 100
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw.text((imWidth / 20, imHeight - imHeight / 20 - fontsize * 1.3), text, fill='white', font=font, anchor=None)
    del draw
    
    if makePicsOfBAYES:
        getProbFuncR = getProbFuncBayes
    else:
        getProbFuncR = getProbFunc
    genes = [i.upper() for i in genes]
    if makePicsOfARGMAX:
        t = Time[np.argmax(np.prod([getProbFuncR(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0))]
    else:
        t = round(np.ceil(np.sum(Time*np.prod([getProbFuncR(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0)/np.sum(np.prod([getProbFunc(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0)))))
    
    text = 'Day '+str(t)
    fontsize = 200
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw = ImageDraw.Draw(img)
    draw.text((imWidth / 20, imHeight / 20), text, fill='white', font=font, anchor=None)
    del draw
    
    os.chdir(dname)
    img.save(fname)
    
    HSV = tuple(map(lambda x: colorsys.hsv_to_rgb(x/len(ims), 1, brightness), range(len(ims))))
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0,0,1,1])
    ax.bar(genes,RNAcount,color=HSV)
    plt.ylabel('Dot Count')
    plt.savefig(fname[:-4]+'_barplot.png', bbox_inches = "tight")
    plt.show()
    
    print(fname+' created!')
    return None

for bmVuc in ['bm','uc']:
    os.chdir(FOLDER7+bmVuc)
    Cells = [f+'/' for f in listdir()]
    Masks = [[c+m for m in listdir(c) if '.tif' in m] for c in Cells]
    Ims = [[c+'gene_nucleus/'+m for m in listdir(c+'gene_nucleus/') if 'MAX' in m][1:] for c in Cells]
    Genes = ['gapdh','col1a1','sox9','runx1','eef2','actb','spp1']
    
    Ims = [[i for i,j in zip(k,Genes) if j.upper() in gene] for k in Ims]
    Genes = list(filter(lambda x: x.upper() in gene, Genes))
    for i,j,k in zip(Cells,Masks,Ims):
        for n,l in enumerate(j):
            os.chdir(FOLDER7+bmVuc)
            makeImage7(l, k, Genes, 'Pics/7_gene/Bar/'+i[:-1]+'_'+str(n)+'.png')





def makeImage(ims, genes, fname):
    if not path.exists(ims[-1][:ims[-1].rfind('/')]+'/'+'cell1_Mask.tif'):
        return None
    rgb = 0
    brightness = 1
    RNAcount = []
    for cc, file in enumerate(ims):
        mask = file[:file.rfind('/')]+'/'+'cell1_Mask.tif'
        M = plt.imread(mask)
        if len(M.shape)==3:
            M=M[:,:,0]
        I = plt.imread(file).astype('float64')
        
        G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
        dd = makePretty(G*(G>0)/np.max(G)*M)
        HSV = colorsys.hsv_to_rgb(cc/len(ims), 1, brightness)
        rgb = rgb + np.concatenate( [dd*HSV[0],dd*HSV[1],dd*HSV[2]],axis=2)
        B = G>cutoff
        RNAcount.append(np.sum(local_maxima(B*G*M)))
    img = Image.fromarray(np.clip(rgb*5,0,255).astype(np.uint8))
    imWidth, imHeight = img.size
    draw = ImageDraw.Draw(img) #2304 pixel = 249.6 micron
    l = (imWidth / 20, imHeight - imHeight / 20,
          imWidth / 20 + round(50*2304/249.6), imHeight - imHeight / 20)
    draw.line(l, fill='white', width=10)
    text = '50 \u03BCm'
    fontsize = 100
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw.text((imWidth / 20, imHeight - imHeight / 20 - fontsize * 1.3), text, fill='white', font=font, anchor=None)
    del draw
    
    
    if makePicsOfBAYES:
        getProbFuncR = getProbFuncBayes
    else:
        getProbFuncR = getProbFunc
    genes = [i.upper() for i in genes]
    if makePicsOfARGMAX:
        t = Time[np.argmax(np.prod([getProbFuncR(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0))]
    else:
        t = round(np.ceil(np.sum(Time*np.prod([getProbFuncR(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0)/np.sum(np.prod([getProbFunc(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0)))))
     
    
    text = 'Day '+str(t)
    fontsize = 200
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw = ImageDraw.Draw(img)
    draw.text((imWidth / 20, imHeight / 20), text, fill='white', font=font, anchor=None)
    del draw
    
    os.chdir(dname)
    img.save(fname)
    
    HSV = tuple(map(lambda x: colorsys.hsv_to_rgb(x/len(ims), 1, brightness), range(len(ims))))
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0,0,1,1])
    ax.bar(genes,RNAcount,color=HSV)
    plt.ylabel('Dot Count')
    plt.savefig(fname[:-4]+'_barplot.png', bbox_inches = "tight")
    plt.show()
    
    print(fname+' created!')
    return None
def getimgs(cell,pathway):
    pathway=pathway.upper()
    gapdh = '1_gapdh_actb/'+cell+'/MAX_C2-'+pathway+'-GAPDH_ACTB_'+cell+'.tif'
    actb = '1_gapdh_actb/'+cell+'/MAX_C4-'+pathway+'-GAPDH_ACTB_'+cell+'.tif'
    il8 = '2_il8_il6_ccl11/'+cell+'/MAX_C2-'+pathway+'-il8_il6_ccl11_'+cell+'.tif'
    il6 = '2_il8_il6_ccl11/'+cell+'/MAX_C3-'+pathway+'-il8_il6_ccl11_'+cell+'.tif'
    ccl11 = '2_il8_il6_ccl11/'+cell+'/MAX_C4-'+pathway+'-il8_il6_ccl11_'+cell+'.tif'
    pathway=pathway.lower()
    col1a1 = '3_col1a1_nanog/'+cell+'/MAX_C2-'+pathway+'-col1a1_nanog_'+cell+'.tif'
    nanog = '3_col1a1_nanog/'+cell+'/MAX_C4-'+pathway+'-col1a1_nanog_'+cell+'.tif'
    sox9 = '4_sox9_eef2_spp1/'+cell+'/MAX_C2-'+pathway+'-sox9_eef2_spp1_'+cell+'.tif'
    eef2 = '4_sox9_eef2_spp1/'+cell+'/MAX_C3-'+pathway+'-sox9_eef2_spp1_'+cell+'.tif'
    spp1 = '4_sox9_eef2_spp1/'+cell+'/MAX_C4-'+pathway+'-sox9_eef2_spp1_'+cell+'.tif'
    runx1 = '5_runx1_pdl1/'+cell+'/MAX_C2-'+pathway+'-runx1_pdl1_'+cell+'.tif'
    pdl1 = '5_runx1_pdl1/'+cell+'/MAX_C4-'+pathway+'-runx1_pdl1_'+cell+'.tif'
    return [gapdh,actb,il8,il6,ccl11,col1a1,nanog,sox9,eef2,spp1,runx1,pdl1]


N = 40
for pathway in ['BM','UC']:
    os.chdir(FOLDER+pathway.lower())
    Cells = [f for f in listdir('5_runx1_pdl1')]
    Ims = [getimgs(c,pathway) for c in Cells]
    Genes = 'gapdh,actb,il8,il6,ccl11,col1a1,nanog,sox9,eef2,spp1,runx1,pdl1'.split(',')
    Ims = [[i for i,j in zip(k,Genes) if j.upper() in gene] for k in Ims]
    Genes = list(filter(lambda x: x.upper() in gene, Genes))
    for i,k in zip(Cells[:N],Ims[:N]):
        os.chdir(FOLDER+pathway.lower())
        makeImage(k, Genes, 'Pics/'+pathway+'/Bar/'+i+'.png')
    

def makeImage3D(mask, ims, genes, fname):
    M = plt.imread(mask)
    rgb = 0
    brightness = 1
    RNAcount = []
    for cc, file in enumerate(ims):
        I = plt.imread(file).astype('float64')
        
        G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
        dd = makePretty(G*(G>0)/np.max(G)*M)
        HSV = colorsys.hsv_to_rgb(cc/len(ims), 1, brightness)
        rgb = rgb + np.concatenate( [dd*HSV[0],dd*HSV[1],dd*HSV[2]],axis=2)
        B = G>cutoff
        RNAcount.append(np.sum(local_maxima(B*G*M)))
    img = Image.fromarray(np.clip(rgb*5,0,255).astype(np.uint8))
    imWidth, imHeight = img.size
    draw = ImageDraw.Draw(img) #2304 pixel = 249.6 micron
    l = (imWidth / 20, imHeight - imHeight / 20,
          imWidth / 20 + round(50*2304/249.6), imHeight - imHeight / 20)
    draw.line(l, fill='white', width=10)
    text = '50 \u03BCm'
    fontsize = 100
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw.text((imWidth / 20, imHeight - imHeight / 20 - fontsize * 1.3), text, fill='white', font=font, anchor=None)
    del draw
    
    
    if makePicsOfBAYES:
        getProbFuncR = getProbFuncBayes
    else:
        getProbFuncR = getProbFunc
    genes = [i.upper() for i in genes]
    if makePicsOfARGMAX:
        t = Time[np.argmax(np.prod([getProbFuncR(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0))]
    else:
        t = round(np.ceil(np.sum(Time*np.prod([getProbFuncR(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0)/np.sum(np.prod([getProbFunc(g,d,RNAcount[genes.index(NormalizationGene)]) for g,d in zip(genes,RNAcount)],axis=0)))))
    
    text = 'Day '+str(t)
    fontsize = 200
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw = ImageDraw.Draw(img)
    draw.text((imWidth / 20, imHeight / 20), text, fill='white', font=font, anchor=None)
    del draw
    
    os.chdir(dname)
    img.save(fname)
    
    
    HSV = tuple(map(lambda x: colorsys.hsv_to_rgb(x/len(ims), 1, brightness), range(len(ims))))
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0,0,1,1])
    ax.bar(genes,RNAcount,color=HSV)
    plt.ylabel('Dot Count')
    plt.savefig(fname[:-4]+'_barplot.png', bbox_inches = "tight")
    plt.show()
    
    print(fname+' created!')
    return None

os.chdir(FOLDER3D)
Cells = [f+'/' for f in listdir() if not '.txt' in f]
Masks = [[c+'count_mask/'+m for m in listdir(c+'count_mask/')] for c in Cells]
Ims = [[c+'nucleus_genes_segmentation/'+m for n,m in enumerate(listdir(c+'nucleus_genes_segmentation/')) if n in (1,5,6,7,8,9,10,11,12,13)] for c in Cells]
Genes = ['pdl1','sox9','cl1a1','actb','il8','il6','ccl11','runx1','spp1','eef2']

Ims = [[i for i,j in zip(k,Genes) if j.upper() in gene] for k in Ims]
Genes = list(filter(lambda x: x.upper() in gene, Genes))

for i,j,k in zip(Cells,Masks,Ims):
    for l in j:
        os.chdir(FOLDER3D)
        makeImage3D(l, k, Genes, 'Pics/3D/Bar/'+i[:-1]+'_'+l[:l.rfind('_')][-1]+'.png')
