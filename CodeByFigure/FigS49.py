#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 15:37:28 2021

@author: ryan
"""


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
import colorsys
import csv
import statistics
import scipy.stats as st
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)



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

"""
coating data
"""
# csv_file = open('JustData.csv')
# read_tsv = csv.reader(csv_file, delimiter=",")
# data = []
# for row in read_tsv:
#     data.append(row)
# csv_file.close()
measure = lambda x: np.log(x)
Data = np.array(list(map(lambda x: np.array(list(map(lambda y: measure(float(y)),x[1:]))), data[23:27])))
fig, ax = plt.subplots()
X = Data[0,:]
Y = Data[2,:]
bp = plt.plot(X,Y,'r.')
X_mean = np.mean(X)
Y_mean = np.mean(Y)
m = sum((X - X_mean)*(Y - Y_mean)) / sum(np.square(X - X_mean))
c = Y_mean - m*X_mean
plt.plot([min(X), max(X)], [min(X)*m+c, max(X)*m+c], 'r--')
stdm = (statistics.variance(Y)/sum(np.square(X - np.mean(X))))**0.5
pb1 = st.norm(m, stdm).cdf(1)
# print('Poly-L-Lysine slope: ' + str(m))
# print('Standard deviation of Poly-L-Lysine slope: ' + str(stdm))
# print('Probability that slope of Poly-L-Lysine is below black line: ' + str(pb1))
X = Data[1,:]
Y = Data[3,:]
bp = plt.plot(X,Y,'b.')
X_mean = np.mean(X)
Y_mean = np.mean(Y)
m = sum((X - X_mean)*(Y - Y_mean)) / sum(np.square(X - X_mean))
c = Y_mean - m*X_mean
plt.plot([min(X), max(X)], [min(X)*m+c, max(X)*m+c], 'b--')
stdm = (statistics.variance(Y)/sum(np.square(X - np.mean(X))))**0.5
pb2 = st.norm(m, stdm).cdf(1)
# print('Collagen slope: ' + str(m))
# print('Standard deviation of Collagen slope: '+str(stdm))
# print('Probability that slope of Collagen is below black line: ' + str(pb2))
plt.plot([min(X), max(X)], [min(X), max(X)], 'k')
ax.set_xlabel(r'$3x10^4/cm^2$ Secretion [log(afu)]')
ax.set_ylabel(r'$10^5/cm^2$ Secretion [log(afu)]')
ax.legend(['Poly-L-Lysine','Poly-L-Lysine Fit','Collagen','Collagen Fit','Independent of Density'])
fig.text(0,0.9,'           Probability that both fit slopes are independently \n           below black line: ' +str(round(pb2*pb1,5)) )
plt.savefig("{}.jpg".format('Quantification/LuminexDensityVSCytokine'))
