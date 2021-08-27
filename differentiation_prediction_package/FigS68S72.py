#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 12:47:13 2021

@author: ryan
"""
import os
from os import path
from os import listdir
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont 
import seaborn as sns
from scipy import ndimage
from skimage.morphology import local_maxima, disk
from functools import reduce
import colorsys
import tensorflow as tf
import random
import tifffile
from scipy.ndimage import correlate
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


"""
HOW THE GENE STUFF IS NORMALIZED FOR ENCLUDING BOTH DATASETS
"""
# NORMIZE = lambda x,y: (np.log(x+0.1)/np.log(y+0.1))
NORMIZE = lambda x,y: np.log( (x+0.1)/(y+0.1) )
# NORMIZE = lambda x,y: (x+0.1)/(y+0.1)
makePicsOfARGMAX = 0
makePicsOfBAYES = 0



# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109503
genesRaw = [i for i in open('PredictionData/ChondrogenesisOnlineData.txt',"r")]

donor = [0,1,2]
Time = np.array([0,1,3,7,14,21])

cats = genesRaw[0].split()
geneind = cats.index('geneName')
gene = [i.split()[geneind] for i in genesRaw[1:]]
timeind = translist([[cats.index(i)-2 for i in cats if 'D'+str(j)+'D' in i] for j in [1,2,3]])

exp = [i.split()[2:] for i in genesRaw[1:]]
exp = tuple(map(lambda  x: [[float(x[k]) for k in l] for l in translist(timeind)], exp))

Exp = {gene[i]: exp[i] for i in range(len(gene))}
    
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
    ProbT = tuple(map(lambda m, s: 0.01 + np.nan_to_num((e**(-0.5*((m - NORMIZE(ValSeen,EValSeen))/s)**2)/s)), mu, sig))
    Porportionality = 1/sum(ProbT)
    return tuple(map(lambda x: x*Porportionality, ProbT))



for cult in ['','3D']:
    with open('RNAcountsBM'+cult+'.txt') as f:
        lines = f.readlines()
    BM = [[int(j) for j in i.split()] for i in lines[1:]]
    tpBM = ['BM' for i in BM]
    with open('RNAcountsUC'+cult+'.txt') as f:
        lines = f.readlines()
    UC = [[int(j) for j in i.split()] for i in lines[1:]]
    tpUC = ['UC' for i in UC]
    
    HCRgene = [i.upper() for i in lines[0].split()]
    
    BM = [ translist(tuple(filter(lambda x: x[0] in gene, zip(HCRgene,i))))[1] for i in BM]
    UC = [ translist(tuple(filter(lambda x: x[0] in gene, zip(HCRgene,i))))[1] for i in UC]
    HCRgene = tuple(filter(lambda x: x in gene, HCRgene))
    
    palette = sns.color_palette()
    
    
    TrainingInfo = translist([translist(Exp[i]) for i in HCRgene])
    TrainingDay = np.array(Time)
    def make1data():
        ThisTimeInd = int(random.random()*len(Time))
        thisRand0 = random.random()
        thisRand1 = random.random()
        thisRand2 = random.random()
        tRR = thisRand0+thisRand1+thisRand2
        Input = [np.random.normal(j[0]*thisRand0/tRR + j[1]*thisRand1/tRR + j[2]*thisRand2/tRR, np.std(j), 1) for j in TrainingInfo[ThisTimeInd]]
        Output = ThisTimeInd
        return [Input, Output]
    train_data = translist([make1data() for i in range(1000*len(HCRgene))])
    x_train = np.array(train_data[0])
    y_train = np.array(train_data[1])
    model = tf.keras.models.Sequential([
      tf.keras.layers.Flatten(input_shape= x_train[0].shape ),
      tf.keras.layers.Dense(20, activation='relu'),
      tf.keras.layers.Dropout(0.2),
      tf.keras.layers.Dense(len(Time), activation='softmax')
    ])
    predictions = model(x_train[:1]).numpy()
    predictions
    tf.nn.softmax(predictions).numpy()
    loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
    loss_fn(y_train[:1], predictions).numpy()
    
    model.compile(optimizer='adam',
                  loss=loss_fn,
                  metrics=['accuracy'])
    model.fit(x_train, y_train, epochs=50)
    
    probability_model = tf.keras.Sequential([
      model,
      tf.keras.layers.Softmax()
    ])
    UCtime = probability_model(np.array(UC))
    BMtime = probability_model(np.array(BM))
    
    EXUC = np.dot(UCtime,Time)
    EXBM = np.dot(BMtime,Time)
    
    bm,x=np.histogram(EXUC,bins=100,range=(min(Time),max(Time)))
    uc,x=np.histogram(EXBM,bins=100,range=(min(Time),max(Time)))
    x=x[:-1]
    ax = plt.subplot(111)
    line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
    line2=ax.bar(x, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
    ax.autoscale(tight=True)
    ax.set_title('Prediction '+cult)
    ax.set_xlabel('Differentiation Time')
    ax.legend([line1, line2], ['BM', 'UC'])
    plt.savefig('Quantification/DifferentiaionPredictionNN'+cult+'.png')
    plt.show()
    
    AMUC = np.argmax(UCtime,axis=1)
    AMBM = np.argmax(BMtime,axis=1)
    
    bm,x=np.histogram(AMUC,bins=100,range=(min(Time),max(Time)))
    uc,x=np.histogram(AMBM,bins=100,range=(min(Time),max(Time)))
    x=x[:-1]
    ax = plt.subplot(111)
    line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
    line2=ax.bar(x, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
    ax.autoscale(tight=True)
    ax.set_title('Prediction '+cult)
    ax.set_xlabel('Differentiation Time')
    ax.legend([line1, line2], ['BM', 'UC'])
    plt.savefig('Quantification/DifferentiaionPredictionNN'+cult+'ARGMAX.png')
    plt.show()
    
    
    BMtime = [ np.sum(Time*np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in BM]
    UCtime = [ np.sum(Time*np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFunc(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in UC]
    bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
    uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
    x=x[:-1]
    ax = plt.subplot(111)
    line1=ax.bar(x, bm/np.sum(bm), width=0.2, alpha=0.5, color=palette[5], align='center')
    line2=ax.bar(x, uc/np.sum(uc), width=0.2, alpha=0.5, color=palette[8], align='center')
    ax.autoscale(tight=True)
    ax.set_title('Prediction '+cult)
    ax.set_xlabel('Differentiation Time')
    ax.legend([line1, line2], ['BM', 'UC'])
    plt.savefig('Quantification/DifferentiaionPredictionSynthetic'+cult+'.png')
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
    ax.set_title('Prediction '+cult)
    ax.set_xlabel('Differentiation Time')
    ax.legend([line1, line2], ['BM', 'UC'])
    plt.savefig('Quantification/DifferentiaionPredictionSynthetic'+cult+'ARGMAX.png')
    plt.show()
    
    BMtime = [ np.sum(Time*np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in BM]
    UCtime = [ np.sum(Time*np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0)/np.sum(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))) for c in UC]
    bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
    uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
    x=x[:-1]
    ax = plt.subplot(111)
    line1=ax.bar(x, bm/np.sum(bm), width=0.2, color=palette[5], align='center')
    line2=ax.bar(x+0.2, uc/np.sum(uc), width=0.2, color=palette[8], align='center')
    ax.autoscale(tight=True)
    ax.set_title('Prediction '+cult)
    ax.set_xlabel('Differentiation Time')
    ax.legend([line1, line2], ['BM', 'UC'])
    plt.savefig('Quantification/DifferentiaionPredictionBayes'+cult+'.png')
    plt.show()
    
    BMtime = [ Time[np.argmax(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))] for c in BM]
    UCtime = [ Time[np.argmax(np.prod([getProbFuncBayes(g,d) for g,d in zip(HCRgene,c)],axis=0))] for c in UC]
    bm,x=np.histogram(BMtime,bins=100,range=(min(Time),max(Time)))
    uc,x=np.histogram(UCtime,bins=100,range=(min(Time),max(Time)))
    x=x[:-1]
    ax = plt.subplot(111)
    line1=ax.bar(x, bm/np.sum(bm), width=0.2, color=palette[5], align='center')
    line2=ax.bar(x+0.2, uc/np.sum(uc), width=0.2, color=palette[8], align='center')
    ax.autoscale(tight=True)
    ax.set_title('Prediction '+cult)
    ax.set_xlabel('Differentiation Time')
    ax.legend([line1, line2], ['BM', 'UC'])
    plt.savefig('Quantification/DifferentiaionPredictionBayes'+cult+'ARGMAX.png')
    plt.show()
    
    print(cult+' Predictions From: '+' '.join(HCRgene))


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



os.chdir('/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/back_sub/001/')
GENES = [f+'/'+j for f in listdir() for j in listdir(f) if (not '_Mask' in j) and (not 'nucleus' in j)]
msks = [j for j in listdir(listdir()[0]) if ('_Mask' in j)]


def JustCount(file,mask):
    I = tifffile.imread(file)
    M = tifffile.imread(file[:file.find('/')+1]+mask)
    gene = file[file.find('/')+1:]
    gene = gene[:gene.find('_')]
    G = ndimage.convolve(I.astype(float),convolvFilter, mode='constant', cval=0)
    # print((np.max(G),np.mean(G),np.min(G)))
    B = G>1000
    DotCount = np.sum(local_maxima(B*G*M))
    return gene, DotCount#, G*(G>0)
# gcs = [[JustCount(file,mask) for file in GENES[-2:-1]] for mask in msks[:1]]
gcs = [[JustCount(file,mask) for file in GENES] for mask in msks]
gcs = [[i for i in j if i[0].upper() in gene] for j in gcs]
gcs = [[i for i in j if i[0].upper() in ('GAPDH','COL1A1','SOX9','RUNX1','EEF2','ACTB','SPP1',)] for j in gcs]
# # I = gcs[0][0][-1]
# # print([[i[:-1] for i in j if i[0].upper() == 'SOX9'] for j in gcs])
# # print(gcs)
# for i in gcs:
#     plt.imshow(binary_dilation(i[0][-1]>1000,disk(20)))
#     plt.show()
# print([i[0] for i in gcs[0] if i[0].upper() in gene])#4
# plt.bar(*translist(gcs[msks.index('cell4_Mask.tif')]))
# plt.xticks(rotation = 45)
# plt.show()
for makePicsOfARGMAX in (0,1):
    for makePicsOfBAYES in (0,1):
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
        print((times, 'bayes', makePicsOfBAYES, 'argmax', makePicsOfARGMAX, ))
# translist(gcs)[-2]

















def makePretty(I):
    I = np.array(I,dtype=float)
    I = (I-np.min(I))/(np.max(I)-np.min(I))*255
    return np.reshape(I, (I.shape[0],I.shape[1],1))

def makeImage3Dfinale(FINALE, fname):
    rgb = 0
    brightness = 1
    geens = []
    cts = []
    for cc, I in enumerate(FINALE):
        ge, DotCount, SAVEPIC = I
        geens.append(ge)
        cts.append(DotCount)
        dd = makePretty(SAVEPIC*(SAVEPIC>0)/np.max(SAVEPIC))
        HSV = colorsys.hsv_to_rgb(cc/len(FINALE), 1, brightness)
        rgb = rgb + np.concatenate( [dd*HSV[0],dd*HSV[1],dd*HSV[2]],axis=2)
        
    if makePicsOfBAYES:
        getProbFuncR = getProbFuncBayes
    else:
        getProbFuncR = getProbFunc
    Pthis = [getProbFuncR(ge,ctss,cts[geens.index(NormalizationGene)]) for ge, ctss in zip(geens,cts)]
    Ptot = np.prod(np.array(Pthis),axis=0)
    img = Image.fromarray(np.clip(rgb*5,0,255).astype(np.uint8))
    imWidth, imHeight = img.size
    draw = ImageDraw.Draw(img) #2304 pixel = 249.6 micron
    l = (imWidth / 160, imHeight - imHeight / 20,
          imWidth / 20 + round(28.17/390 *1000), imHeight - imHeight / 20)
    draw.line(l, fill='white', width=10)
    text = '1000 \u03BCm'
    fontsize = 10
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw.text((imWidth / 20, imHeight - imHeight / 20 - fontsize * 1.3), text, fill='white', font=font, anchor=None)
    del draw
    if makePicsOfARGMAX:
        t = Time[np.argmax(Ptot)]
    else:
        t = round(np.ceil(np.sum(Time*Ptot/np.sum(Ptot))))
    text = 'Day '+str(t)
    fontsize = 20
    font = ImageFont.truetype('Arial.ttf', fontsize)
    draw = ImageDraw.Draw(img)
    draw.text((imWidth / 20, imHeight / 20), text, fill='white', font=font, anchor=None)
    del draw
    os.chdir(dname)
    img.save(fname)
    HSV = tuple(map(lambda x: colorsys.hsv_to_rgb(x/len(FINALE), 1, brightness), range(len(FINALE))))
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0,0,1,1])
    ax.bar(geens,cts,color=HSV)
    plt.ylabel('Dot Count')
    plt.savefig(fname[:-4]+'_barplot.png', bbox_inches = "tight")
    plt.show()
    print(fname+' created!')
    return None
os.chdir('/media/ryan/Seagate Expansion Drive/DeepLearningRNAscope/back_sub/')
GENES = [f+'/' for f in listdir()]
geness = [(int(file[file.find('-')-1]),file[:file.find('_')].upper(),file[:-1] )for file in GENES if (file[:file.find('_')].upper() in gene) and (not file[:file.find('_')]=='nucleus') and ('.tif' in file)]

def ReturnDotsAndPic(totinfo):
    ch, ge, fi = totinfo
    ims = tifffile.imread(fi)
    SAVEPIC = np.max(ims, axis=0)
    G = ndimage.convolve(SAVEPIC,convolvFilter, mode='constant', cval=0)
    n1 = 12
    n2 = 11
    n3 = 6
    n4 = 5
    h1 = 9
    h2 = 3
    convolvFilter0 = -disk(n1).astype(float)
    convolvFilter0[n1-n2:len(convolvFilter0)-n1+n2,n1-n2:len(convolvFilter0)-n1+n2] += disk(n2).astype(float)
    convolvFilter0[n1-n3:len(convolvFilter0)-n1+n3,n1-n3:len(convolvFilter0)-n1+n3] -= disk(n3).astype(float)*np.sum(convolvFilter0)/np.sum(disk(n3))
    convolvFilter1 = np.array(convolvFilter0)
    convolvFilter2 = np.array(convolvFilter0)
    convolvFilter2[n1-n4:len(convolvFilter0)-n1+n4,n1-n4:len(convolvFilter0)-n1+n4] += disk(n4).astype(float)*np.sum(disk(n1))/np.sum(disk(n4))*2/h2
    t = [-disk(n1).astype(float)[np.newaxis,:,:]]
    m = int((h1-h2)/2)*[convolvFilter1[np.newaxis,:,:]]
    c = h2*[convolvFilter2[np.newaxis,:,:]]
    convolvFilter3D = np.vstack(t+m+c+m+t)
    convolvFilter3D = convolvFilter3D/np.sum(abs(convolvFilter3D))
    
    DotProb = np.array(correlate(ims[:,:,:],convolvFilter3D))
    DotLoc = np.array(local_maxima((DotProb > 50000)*ims) )
    DotCount = np.sum(DotLoc)
    return [ge, DotCount, SAVEPIC]
FINALE = [ReturnDotsAndPic(i) for i in geness]
FinaleSample = ['col1a1'.upper(), 106,
                'eef2'.upper(), 804,
                'spp1'.upper(), 68,#1/4 to 1/3 of sox9 or less
                'il8'.upper(), 62,
                'il6'.upper(), 48,
                'ccl11'.upper(), 37,
                'runx1'.upper(), 156,
                'pdl1'.upper(), 3,
                'sox9'.upper(), 62,
                'actb'.upper(), 45,]
FINALE = [ [i[0],  FinaleSample[1::2][FinaleSample[::2].index(i[0])]]  +  [i[2]] for i in FINALE if i[0] in FinaleSample]
makeImage3Dfinale(FINALE,'Pics/finale/back_sub.png')





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

