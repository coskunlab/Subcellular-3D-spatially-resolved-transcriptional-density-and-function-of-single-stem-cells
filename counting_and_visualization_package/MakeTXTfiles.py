#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:34:02 2021

@author: ryan
"""

import os
from os import path
from os import listdir
from os.path import isdir
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from skimage.morphology import local_maxima, disk, binary_dilation, binary_erosion
import itertools
import pandas as pd

dname = os.path.dirname(os.path.abspath(__file__)) #PYTHON SCRIPT DIRECTORY


"""
COUNTING SECTION 05-12

THIS GENERATES .txt FILES OF THE LOCATIONS OF ALL RNAs IN THE 05-12 DATASET
"""
FOLDER = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/12_gene_data_05-12/'
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
def flatten(t):
    return [item for sublist in t for item in sublist]

def CountDots(mask,file):
    M = plt.imread(mask)
    if len(M.shape) == 3:
        M = M[:,:,0]
    I = plt.imread(file).astype('float64')
    G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
    B = (G>cutoff)
    D = local_maxima(B*G*M)
    return np.sum(D)

Cells = [(3-len(str(i)))*'0'+str(i)+'/' for i in range(1,46)]
for typ in [FOLDER+'bm/',FOLDER+'uc/']:
    for dns in [typ+'high_dense/',typ+'low_dense/']:
        if dns == FOLDER+'uc/'+'high_dense/':
            Get_Well = lambda x: 1 + 1*(x in range(24,36))
        else:
            Get_Well = lambda x: np.ceil(x/15)
        os.chdir(FOLDER)
        pos = pd.read_excel(FOLDER+'positions.xlsx', sheet_name= (typ[76:78].upper()+'_'+dns[79:][:-1].capitalize()) )
        os.chdir(dns)
        Masks = [[c+m for m in listdir(c) if 'Mask' in m] for c in Cells if isdir(c)]
        Ims = [[c+m for m in listdir(c) if 'MAX' in m and not 'C1' in m] for c in Cells if isdir(c)]
        CellNum = [i[0][:3] for i in Ims]
        Cols = ['well','x','y','il8','il6','ccl11']
        alldotcounts=flatten([[ [Get_Well(int(i)),pos.x[int(i)-1],pos.y[int(i)-1]]+[CountDots(m,l) for l in k] for m in j] for i,j,k in zip(CellNum,Masks,Ims)])
        os.chdir(dname)
        file1 = open("RNAcounts05_12"+typ[76:78].upper()+dns[79:][:-1]+".txt","w")
        GG = ' '.join(Cols)
        CC = '\n'.join([' '.join([str(j) for j in i]) for i in alldotcounts])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close() 
        print('Did a thing!')




"""
COUNTING SECTION 3D

THIS GENERATES .txt FILES OF THE COUNTS OF EACH CELL'S RNAs IN THE 3D DATASET
"""
FOLDER3D = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/3D dataset/'
FOLDER3Dpt2 = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/give Ryan/'

def translist(x):
    return [list(i) for i in zip(*x)]
def ListMult(x,y):
    return [[i,j] for i in x for j in y]


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
def flatten(t):
    return [item for sublist in t for item in sublist]

def CountDots3D(mask,file):
    M = plt.imread(mask)
    if len(M.shape) == 3:
        M = M[:,:,0]
    I = plt.imread(file).astype('float64')
    G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
    B = (G>cutoff)
    D = local_maxima(B*G*M)
    if ('il6_MAX' in file) and ('3D dataset' in file):
        plt.imshow(I,cmap='gray')
        plt.show()
        plt.imshow(G*(G>0),cmap='gray')
        plt.show()
        plt.imshow(binary_dilation(D,disk(20)))
        plt.show()
        plt.imshow(M)
        plt.show()
        print(file)
        print(np.sum(D))
    return np.sum(D)
os.chdir(FOLDER3D)
Cells = [f+'/' for f in listdir() if not '.txt' in f]
Masks = [[c+'count_mask/'+m for m in listdir(c+'count_mask/')] for c in Cells]
Ims = [[c+'nucleus_genes_segmentation/'+m for n,m in enumerate(listdir(c+'nucleus_genes_segmentation/')) if n in (1,5,6,7,8,9,10,11,12,13)] for c in Cells]
Ims = [sorted(i, key=lambda x: x[37:40]) for i in Ims]
Genes = ['ccl11','actb','eef2','il6','il8','col1a1','sox9','pdl1','spp1','runx1']
alldotcounts=flatten([[[k[0][:2]]+[CountDots3D(m,l) for l in k] for m in j] for j,k in zip(Masks,Ims)])
BM = [i[1:] for i in alldotcounts if i[0]=='bm']
UC = [i[1:] for i in alldotcounts if i[0]=='uc']

os.chdir(FOLDER3Dpt2)
Cells = [f+'/' for f in listdir()]
MIG= [i for i in flatten([[flatten([ListMult([c+g+'/'+m for m in listdir(c+'/'+g) if 'Mask' in m and str(mnum) in m],[(c+g+'/'+m,m[:m.find('_')]) for m in listdir(c+g) if not 'Mask' in m]) for g in listdir(c)]) for c in Cells] for mnum in range(5)]) if i]
Counts = [[(CountDots3D(f[0],f[1][0]),f[1][1]) for f in c] for c in MIG]
UCnew = [translist(sorted(c, key=lambda x: Genes.index(x[1])))[0] for c in Counts]

UC = UC + UCnew

os.chdir(dname)
file1 = open("RNAcountsBM3D.txt","w")
GG = ' '.join(Genes)
CC = '\n'.join([' '.join([str(j) for j in i]) for i in BM])
file1.write(GG+'\n')
file1.writelines(CC)
file1.close() 
file1 = open("RNAcountsUC3D.txt","w")
GG = ' '.join(Genes)
CC = '\n'.join([' '.join([str(j) for j in i]) for i in UC])
file1.write(GG+'\n')
file1.writelines(CC)
file1.close() 




"""
COUNTING SECTION

THIS GENERATES .txt FILES OF THE COUNTS OF EACH CELL'S RNAs IN THE 2D DATASET
"""
FOLDER = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/12_gene_data/'

def translist(list_here):
    return [list(x) for x in zip(*list_here)]

def flatten(t):
    return [item for sublist in t for item in sublist]

os.chdir(FOLDER)

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

Cells = [str(i+1) for i in range(95)]
Cells = ['0'*(3-len(i))+i+'/' for i in Cells]

def CountDots(file,cell):
    M = plt.imread(file[:file.rfind('/')]+'/cell'+str(cell)+'_Mask.tif')
    if len(M.shape) == 3:
        M = M[:,:,0]
    I = plt.imread(file).astype('float64')
    G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
    B = (G>cutoff)
    D = local_maxima(B*G*M)
    return np.sum(D)
def CellCounts(subFOLDER,cell):
    pathway = subFOLDER[-3:]
    gapdh = subFOLDER+'1_gapdh_actb/'+cell+'MAX_C2-'+pathway.upper()[:-1]+'-GAPDH_ACTB_'+str(cell)[:-1]+'.tif'
    actb = subFOLDER+'1_gapdh_actb/'+cell+'MAX_C4-'+pathway.upper()[:-1]+'-GAPDH_ACTB_'+str(cell)[:-1]+'.tif'
    
    il8 = subFOLDER+'2_il8_il6_ccl11/'+cell+'MAX_C2-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    il6 = subFOLDER+'2_il8_il6_ccl11/'+cell+'MAX_C3-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    ccl11 = subFOLDER+'2_il8_il6_ccl11/'+cell+'MAX_C4-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    
    col1a1 = subFOLDER+'3_col1a1_nanog/'+cell+'MAX_C2-'+pathway[:-1]+'-col1a1_nanog_'+str(cell)[:-1]+'.tif'
    nanog = subFOLDER+'3_col1a1_nanog/'+cell+'MAX_C4-'+pathway[:-1]+'-col1a1_nanog_'+str(cell)[:-1]+'.tif'

    sox9 = subFOLDER+'4_sox9_eef2_spp1/'+cell+'MAX_C2-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    eef2 = subFOLDER+'4_sox9_eef2_spp1/'+cell+'MAX_C3-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    spp1 = subFOLDER+'4_sox9_eef2_spp1/'+cell+'MAX_C4-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    
    runx1 = subFOLDER+'5_runx1_pdl1/'+cell+'MAX_C2-'+pathway[:-1]+'-runx1_pdl1_'+str(cell)[:-1]+'.tif'
    pdl1 = subFOLDER+'5_runx1_pdl1/'+cell+'MAX_C4-'+pathway[:-1]+'-runx1_pdl1_'+str(cell)[:-1]+'.tif'
    
    num = 1
    mask = subFOLDER+'5_runx1_pdl1/'+cell+'cell'+str(num)+'_Mask.tif'
    while path.exists(mask):
        num = num+1
        mask = mask[:-10]+str(num)+mask[-9:]
    num = num-1
    return [[CountDots(i,c+1) for i in [gapdh,actb,il8,il6,ccl11,col1a1,nanog,sox9,eef2,spp1,runx1,pdl1]] for c in range(num)]
for pathway in ['bm/','uc/']:
    subFOLDER = FOLDER + pathway
    os.chdir(subFOLDER)
    Dots = [CellCounts(subFOLDER,i) for i in Cells]
    Dots = list(itertools.chain(*Dots))
    Genes = ['gapdh','actb','il8','il6','ccl11','col1a1','nanog','sox9','eef2','spp1','runx1','pdl1']
    os.chdir(dname)
    if pathway=='bm/':
        file1 = open("RNAcountsBM.txt","w")
        GG = ' '.join(Genes)
        CC = '\n'.join([' '.join([str(j) for j in i]) for i in Dots])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close() 
    if pathway=='uc/':
        file1 = open("RNAcountsUC.txt","w")
        GG = ' '.join(Genes)
        CC = '\n'.join([' '.join([str(j) for j in i]) for i in Dots])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close() 


def GetLocs(file,cell):
    M = plt.imread(file[:file.rfind('/')]+'/cell'+str(cell)+'_Mask.tif')
    if len(M.shape) == 3:
        M = M[:,:,0]
    I = plt.imread(file).astype('float64')
    G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
    B = (G>cutoff)#*(I>50) #1e5
    D = local_maxima(B*G*M)
    
    L = np.where(D) #D[P[0],P[1]]
    C = np.mean(np.where(M),axis=1).reshape((2,1)) #D[P[0],P[1]]
    Dc = np.mean(np.sum((L-C)**2,axis=0)**0.5)
    
    Pm = 0
    cou = 0
    while not (np.max(M) == 0):
        cou += 1
        M = binary_erosion(M,disk(10))
        Pm = Pm + M
    dist = np.array((Pm+0.5)*D)
    De = np.sum(dist)*10/np.sum(D)
    return (Dc,De)
def CellAvgLocs(subFOLDER,cell):
    pathway = subFOLDER[-3:]
    gapdh = subFOLDER+'1_gapdh_actb/'+cell+'MAX_C2-'+pathway.upper()[:-1]+'-GAPDH_ACTB_'+str(cell)[:-1]+'.tif'
    actb = subFOLDER+'1_gapdh_actb/'+cell+'MAX_C4-'+pathway.upper()[:-1]+'-GAPDH_ACTB_'+str(cell)[:-1]+'.tif'
    
    il8 = subFOLDER+'2_il8_il6_ccl11/'+cell+'MAX_C2-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    il6 = subFOLDER+'2_il8_il6_ccl11/'+cell+'MAX_C3-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    ccl11 = subFOLDER+'2_il8_il6_ccl11/'+cell+'MAX_C4-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    
    col1a1 = subFOLDER+'3_col1a1_nanog/'+cell+'MAX_C2-'+pathway[:-1]+'-col1a1_nanog_'+str(cell)[:-1]+'.tif'
    nanog = subFOLDER+'3_col1a1_nanog/'+cell+'MAX_C4-'+pathway[:-1]+'-col1a1_nanog_'+str(cell)[:-1]+'.tif'

    sox9 = subFOLDER+'4_sox9_eef2_spp1/'+cell+'MAX_C2-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    eef2 = subFOLDER+'4_sox9_eef2_spp1/'+cell+'MAX_C3-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    spp1 = subFOLDER+'4_sox9_eef2_spp1/'+cell+'MAX_C4-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    
    runx1 = subFOLDER+'5_runx1_pdl1/'+cell+'MAX_C2-'+pathway[:-1]+'-runx1_pdl1_'+str(cell)[:-1]+'.tif'
    pdl1 = subFOLDER+'5_runx1_pdl1/'+cell+'MAX_C4-'+pathway[:-1]+'-runx1_pdl1_'+str(cell)[:-1]+'.tif'
    print('Running: '+subFOLDER[1:]+cell[:-1])
    num = 1
    mask = subFOLDER+'5_runx1_pdl1/'+cell+'cell'+str(num)+'_Mask.tif'
    while path.exists(mask):
        num = num+1
        mask = mask[:-10]+str(num)+mask[-9:]
    num = num-1
    return [[GetLocs(i,c+1) for i in [gapdh,actb,il8,il6,ccl11,col1a1,nanog,sox9,eef2,spp1,runx1,pdl1]] for c in range(num)]

for pathway in ['bm/','uc/']:
    subFOLDER = FOLDER + pathway
    os.chdir(subFOLDER)
    Dots = [CellAvgLocs(subFOLDER,i) for i in Cells]
    Dots = list(itertools.chain(*Dots))
    Genes = ['gapdh','actb','il8','il6','ccl11','col1a1','nanog','sox9','eef2','spp1','runx1','pdl1']
    os.chdir(dname)
    if pathway=='bm/':
        file1 = open("RNAcdistBM.txt","w")
        GG = ' '.join(Genes)
        CC = '\n'.join([' '.join([str(j[0]) for j in i]) for i in Dots])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close() 
        file1 = open("RNAedistBM.txt","w")
        GG = ' '.join(Genes)
        CC = '\n'.join([' '.join([str(j[1]) for j in i]) for i in Dots])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close()
    if pathway=='uc/':
        file1 = open("RNAcdistUC.txt","w")
        GG = ' '.join(Genes)
        CC = '\n'.join([' '.join([str(j[0]) for j in i]) for i in Dots])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close() 
        file1 = open("RNAedistUC.txt","w")
        GG = ' '.join(Genes)
        CC = '\n'.join([' '.join([str(j[1]) for j in i]) for i in Dots])
        file1.write(GG+'\n')
        file1.writelines(CC)
        file1.close() 

