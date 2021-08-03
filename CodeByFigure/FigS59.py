#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 16:06:27 2021

@author: ryan
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib


font_size_of_the_code = 12#24
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : font_size_of_the_code}
matplotlib.rc('font', **font)






"""
3D EXPRESSION AND LOCATION PLOTS
"""
for typ in ['BM','UC']:
    for den in ['high','low']:
        file = 'RNAcounts05_12'+typ+den+'_dense.txt'
        with open(file) as f:
            lines = f.readlines()
        cd = np.array([[float(j) for j in i.split()] for i in lines[1:]]).T
        gene = lines[0].split()
        w = cd[0]
        x = cd[1]/1000
        y = cd[2]/1000
        il8 = cd[3]
        il6 = cd[4]
        ccl11 = cd[5]
        ax = plt.axes(projection='3d')
        ax.scatter3D(x, y, il8, c='yellow')
        ax.scatter3D(x, y, il6, c='cyan')
        ax.scatter3D(x, y, ccl11, c='magenta')
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('y [mm]')
        ax.set_zlabel('Dots')
        ax.legend(gene[3:])
        # ax.plot_trisurf(x, y, il6, alpha=0.5, edgecolor='black')
        # plt.show()
        # ax.scatter3D(x, y, il8, c=il8, cmap='plasma');
        # ax.plot_trisurf(x, y, il8, cmap='magma', edgecolor='none');
        ax.view_init(30, 0)
        plt.savefig('Quantification/location_plot_'+file[14:][:-4])
        plt.show()
        
        a = 0.8
        for well in set(w):
            x = [i for i,j in zip(cd[1]/1000,w) if j == well]
            y = [i for i,j in zip(cd[2]/1000,w) if j == well]
            il8 = [(i) for i,j in zip(cd[3],w) if j == well]
            il6 = [(i) for i,j in zip(cd[4],w) if j == well]
            ccl11 = [(i) for i,j in zip(cd[5],w) if j == well]
            fig, axs = plt.subplots(2,2)
            axs[0,0].title.set_text('Well '+str(well)[0]+': '+file[14:][:-4])
            axs[0,0].scatter(x, y, c=il8, cmap='inferno',alpha=a)
            axs[0,1].scatter(x, y, c=il6, cmap='inferno',alpha=a)
            axs[1,0].scatter(x, y, c=ccl11, cmap='inferno',alpha=a)
            axs[1,1].scatter(np.linspace(np.min(x),np.max(x),100)*0, np.linspace(0,1,100), c=np.linspace(np.min(x),np.max(x),100), cmap='inferno',alpha=a)
            axs[1,0].set_xlabel('location [mm]')
            axs[0,0].set_ylabel('il8')
            axs[0,1].set_ylabel('il6')
            axs[1,0].set_ylabel('ccl11')
            axs[1,1].set_xlabel('colorkey')
            axs[1,1].set_ylabel('dots')
            axs[1,1].get_xaxis().set_ticks([])
            # axs[1,1].get_yaxis().set_ticks([])
            # axs[0,0].set_ylim([-20, 20])
            # axs[0,1].set_ylim([-20, 20])
            # axs[1,0].set_ylim([-20, 20])
            # axs[0,0].set_xlim([np.mean(x)-20, np.mean(x)+20])
            plt.tight_layout()
            plt.savefig('Quantification/location_plot_2D_well'+str(well)+'_'+file[14:][:-4]+'.png')
            plt.show()
            
            
            
            ax = plt.axes(projection='3d')
            ax.plot_trisurf(x, y, il8, alpha=0.5, color='yellow')
            ax.set_xlabel('x [mm]')
            ax.set_ylabel('y [mm]')
            ax.set_zlabel('Dots')
            ax.view_init(30, 45)
            plt.savefig('Quantification/location_well'+str(well)+'_plot_'+file[14:][:-4]+'_il8.png')
            plt.show()
            
            ax = plt.axes(projection='3d')
            ax.plot_trisurf(x, y, il6, alpha=0.5, color='cyan')
            ax.set_xlabel('x [mm]')
            ax.set_ylabel('y [mm]')
            ax.set_zlabel('Dots')
            ax.view_init(30, 45)
            plt.savefig('Quantification/location_well'+str(well)+'_plot_'+file[14:][:-4]+'_il6.png')
            plt.show()
            
            ax = plt.axes(projection='3d')
            ax.plot_trisurf(x, y, ccl11, alpha=0.5, color='magenta')
            ax.set_xlabel('x [mm]')
            ax.set_ylabel('y [mm]')
            ax.set_zlabel('Dots')
            ax.view_init(30, 45)
            plt.savefig('Quantification/location_well'+str(well)+'_plot_'+file[14:][:-4]+'_ccl11.png')
            plt.show()
            
            ax = plt.axes(projection='3d')
            ax.scatter3D(x, y, il8, alpha=0.5, color='yellow')
            ax.set_xlabel('x [mm]')
            ax.set_ylabel('y [mm]')
            ax.set_zlabel('Dots')
            ax.view_init(30, 45)
            plt.savefig('Quantification/location_well'+str(well)+'_plot_'+file[14:][:-4]+'_il8Scatter.png')
            plt.show()
            
            ax = plt.axes(projection='3d')
            ax.scatter3D(x, y, il6, alpha=0.5, color='cyan')
            ax.set_xlabel('x [mm]')
            ax.set_ylabel('y [mm]')
            ax.set_zlabel('Dots')
            ax.view_init(30, 45)
            plt.savefig('Quantification/location_well'+str(well)+'_plot_'+file[14:][:-4]+'_il6Scatter.png')
            plt.show()
            
            ax = plt.axes(projection='3d')
            ax.scatter3D(x, y, ccl11, alpha=0.5, color='magenta')
            ax.set_xlabel('x [mm]')
            ax.set_ylabel('y [mm]')
            ax.set_zlabel('Dots')
            ax.view_init(30, 45)
            plt.savefig('Quantification/location_well'+str(well)+'_plot_'+file[14:][:-4]+'_ccl11Scatter.png')
            plt.show()
            
        