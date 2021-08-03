#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:34:52 2021

@author: ryan
"""

import os
from os import path
from os import listdir
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont 
from scipy import ndimage
import colorsys
from skimage.morphology import local_maxima, disk, binary_dilation
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)


dname = os.path.dirname(os.path.abspath(__file__)) #PYTHON SCRIPT DIRECTORY

FOLDER = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/12_gene_data/'
FOLDER3D = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/3D dataset/'
FOLDER3Dfinale = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/back_sub/'
FOLDERFIG = '/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/cytokine figure images/'
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
    if np.max(I) > 0:
        I = (I-np.min(I))/(np.max(I)-np.min(I))*255
    return np.reshape(I, (I.shape[0],I.shape[1],1))


def makeImageFIG(mask, ims, fname, ds=10):
    M = plt.imread(mask)
    rgb = 0
    brightness = 1
    for cc, file in enumerate(ims):
        I = plt.imread(file).astype('float64')
        # plt.imshow(I)
        # plt.show()
        G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
        B = G>cutoff
        D = local_maxima(B*G*M)
        dd = makePretty(binary_dilation(D,disk(ds)))
        HSV = colorsys.hsv_to_rgb(cc/len(ims)+0.5/len(ims), 1, brightness)
        rgb = rgb + np.concatenate( [dd*HSV[0],dd*HSV[1],dd*HSV[2]],axis=2)
    # plt.imshow(rgb)
    # plt.show()
    
    M = makePretty(M)
    rgb = rgb + M/4
    img = Image.fromarray(np.clip(rgb,0,255).astype(np.uint8))
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
    
    os.chdir(dname)
    img.save(fname)
    print(fname+' created!')
    return None
os.chdir(FOLDERFIG)
Picccs = [f+'/' for f in listdir() if not ('.url' in f)]
Ims = [[p+i for i in listdir(p) if ('MAX' in i) and not ('C1' in i)] for p in Picccs][3:]
Mask = [[p+i for i in listdir(p) if not ('max' in i.lower())][0] for p in Picccs][3:]
for j,i in zip(Mask,Ims):
    os.chdir(FOLDERFIG)
    makeImageFIG(j,i, 'Pics/FIG_PICS/'+i[0][:i[0].find('/')]+'.png')

def makeRGB(fpath,fname,c,y,m=None,ds=10):
    if m == None:
        rgb = 0
        for cc, file in enumerate([c,y]): #r,m
            M = plt.imread(file[:file.rfind('/')]+'/cell1_Mask.tif')
            if len(M.shape) == 3:
                M = M[:,:,0]
            I = plt.imread(file).astype('float64')
            G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
            B = G>cutoff
            D = local_maxima(B*G*M)
            
            mm = makePretty(M)
            dd = makePretty(binary_dilation(D,disk(ds)))
            rgb = rgb + np.concatenate( [dd*(cc!=0),dd*(cc!=1),dd*(cc!=0)],axis=2)+mm/12
    else:
        rgb = 0
        for cc, file in enumerate([c,y,m]): #c,y,m
            M = plt.imread(file[:file.rfind('/')]+'/cell1_Mask.tif')
            if len(M.shape) == 3:
                M = M[:,:,0]
            I = plt.imread(file).astype('float64')
            
            G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
            B = G>cutoff
            D = local_maxima(B*G*M)
            
            mm = makePretty(M)
            dd = makePretty(binary_dilation(D,disk(ds)))
            rgb = rgb + np.concatenate( [dd*(cc!=1),dd*(cc!=2),dd*(cc!=0)],axis=2)+mm/12
    img = Image.fromarray(np.clip(rgb,0,255).astype(np.uint8))
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
    os.chdir(fpath)
    img.save(fname)
    return None
def makeImage(ctype='UC',cell=79):
    fpath=dname+'/Pics'+'/'+ctype.upper()
    if ctype.upper() == 'UC':
        pathway='uc/'
    elif ctype.upper() == 'BM':
        pathway='bm/'
    else:
        print("We don't have that cell type")
        return None
    cell = '0'*(3-len(str(cell)))+str(cell)+'/'
    gapdh = '1_gapdh_actb/'+cell+'MAX_C2-'+pathway.upper()[:-1]+'-GAPDH_ACTB_'+str(cell)[:-1]+'.tif'
    actb = '1_gapdh_actb/'+cell+'MAX_C4-'+pathway.upper()[:-1]+'-GAPDH_ACTB_'+str(cell)[:-1]+'.tif'
    il8 = '2_il8_il6_ccl11/'+cell+'MAX_C2-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    il6 = '2_il8_il6_ccl11/'+cell+'MAX_C3-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    ccl11 = '2_il8_il6_ccl11/'+cell+'MAX_C4-'+pathway.upper()[:-1]+'-il8_il6_ccl11_'+str(cell)[:-1]+'.tif'
    col1a1 = '3_col1a1_nanog/'+cell+'MAX_C2-'+pathway[:-1]+'-col1a1_nanog_'+str(cell)[:-1]+'.tif'
    nanog = '3_col1a1_nanog/'+cell+'MAX_C4-'+pathway[:-1]+'-col1a1_nanog_'+str(cell)[:-1]+'.tif'
    sox9 = '4_sox9_eef2_spp1/'+cell+'MAX_C2-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    eef2 = '4_sox9_eef2_spp1/'+cell+'MAX_C3-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    spp1 = '4_sox9_eef2_spp1/'+cell+'MAX_C4-'+pathway[:-1]+'-sox9_eef2_spp1_'+str(cell)[:-1]+'.tif'
    runx1 = '5_runx1_pdl1/'+cell+'MAX_C2-'+pathway[:-1]+'-runx1_pdl1_'+str(cell)[:-1]+'.tif'
    pdl1 = '5_runx1_pdl1/'+cell+'MAX_C4-'+pathway[:-1]+'-runx1_pdl1_'+str(cell)[:-1]+'.tif'
    
    
    os.chdir(FOLDER+pathway)
    if not path.exists('5_runx1_pdl1/'+cell+'cell1_Mask.tif'):
        print('masks dont exist!')
        return None
    else:
        os.chdir(FOLDER+pathway)
        makeRGB(fpath,'RNAdots'+ctype.upper()+cell[:-1]+'_gapdhactb.png',gapdh,actb,ds=3)
        
        os.chdir(FOLDER+pathway)
        makeRGB(fpath,'RNAdots'+ctype.upper()+cell[:-1]+'_il6il8ccl11.png',il8,il6,ccl11)
        
        os.chdir(FOLDER+pathway)
        makeRGB(fpath,'RNAdots'+ctype.upper()+cell[:-1]+'_col1a1nanog.png',col1a1,nanog)
        
        os.chdir(FOLDER+pathway)
        makeRGB(fpath,'RNAdots'+ctype.upper()+cell[:-1]+'_sox9eef2spp1.png',sox9,eef2,spp1)
        
        os.chdir(FOLDER+pathway)
        makeRGB(fpath,'RNAdots'+ctype.upper()+cell[:-1]+'_runx1pdl1.png',runx1,pdl1)
        
        print(ctype+str(cell)+' pics made!')
        return None
N = 20
for ctype in ['uc','bm']:
    for cell in range(N):
        pass
        makeImage(ctype=ctype,cell=cell)

def makeImage3D(mask, ims, fname, ds=10):
    M = plt.imread(mask)
    rgb = 0
    brightness = 1
    for cc, file in enumerate(ims):
        I = plt.imread(file).astype('float64')
        
        G = ndimage.convolve(I,convolvFilter, mode='constant', cval=0)
        B = G>cutoff
        D = local_maxima(B*G*M)
        dd = makePretty(binary_dilation(D,disk(ds)))
        HSV = colorsys.hsv_to_rgb(cc/len(ims), 1, brightness)
        rgb = rgb + np.concatenate( [dd*HSV[0],dd*HSV[1],dd*HSV[2]],axis=2)
    M = makePretty(M)
    rgb = rgb + M/4
    img = Image.fromarray(np.clip(rgb,0,255).astype(np.uint8))
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
    
    os.chdir(dname)
    img.save(fname)
    print(fname+' created!')
    return None
os.chdir(FOLDER3D)
Cells = [f+'/' for f in listdir() if not '.txt' in f]
Masks = [[c+'count_mask/'+m for m in listdir(c+'count_mask/')] for c in Cells]
Ims = [[c+'nucleus_genes_segmentation/'+m for n,m in enumerate(listdir(c+'nucleus_genes_segmentation/')) if n in (1,5,6,7,8,9,10,11,12,13)] for c in Cells]
Ims = [sorted(i, key=lambda x: x[37:40]) for i in Ims]
Genes = ['ccl11','actb','eef2','il6','il8','col1a1','sox9','pdl1','spp1','runx1']
tp = [i[0][:2] for i in Ims]
for i,j,k in zip(Cells,Masks,Ims):
    for l in j:
        os.chdir(FOLDER3D)
        makeImage3D(l, k, 'Pics/3D/'+i[:-1]+'_'+l[:l.rfind('_')][-1]+'.png', ds=2)