#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 14:51:50 2021

@author: ryan
"""
import matplotlib.pyplot as plt
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)






def translist(x):
    return [list(i) for i in zip(*x)]

colors = ['skyblue','navy','yellow','cyan','magenta','green','lightpink','purple','gold','indigo','violet','crimson']
for ctype in ('BM','UC'):
    with open('RNAedist'+ctype+'.txt') as f:
        lines = f.readlines()
    ed = translist([[float(j) for j in i.split()] for i in lines[1:]])
    with open('RNAcdist'+ctype+'.txt') as f:
        lines = f.readlines()
    cd = translist([[float(j) for j in i.split()] for i in lines[1:]])
    gene = lines[0].split()
    
    fig, ax = plt.subplots(figsize=(8, 10))
    for gn,g in enumerate(gene):
        if g in ('ccl11'):
            s,a = 50,0.9
        elif g in ('il8','il6'):
            s,a = 20,0.9
        else:
            s,a = 20,0.1
        plt.scatter(ed[gn],cd[gn],alpha=a,s=s,c=colors[gn])
    plt.legend(gene)
    plt.xlabel('Edge Distance')
    plt.ylabel('Center Distance')
    plt.title(ctype)
    plt.savefig('Quantification/LocationScatterplot'+ctype+'.png')
    plt.show()