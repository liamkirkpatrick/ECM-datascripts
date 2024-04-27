#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:39:20 2024

Plot all data from ALHIC2302 Deep Cores

@author: Liam
"""

#%% Import packages 

# general
import numpy as np
import pandas as pd

# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

#%% User Inputs

# smoothing window, mm
window = 10

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/basic_plot_2302/'
metadata_file = 'metadata.csv'

# ONE BIG PLOT
onebigplot = True
indplots = True

#%% Read in metadata and import data

meta = pd.read_csv(path_to_data+metadata_file)

# import each script as an ECM class item
data = []
cores = []
sections = []
faces = []
ACorDCs = []
for index,row in meta.iterrows():
    
    core = row['core']
    
    # filter for ALHIC2302 shallow ice
    if core == 'alhic2302' and row['idx_abs'] > 0:

        
        section = row['section']
        face = row['face']
        ACorDC = row['ACorDC']
        
        print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
    
        data_item = ECM(core,section,face,ACorDC)
        data_item.rem_ends(10)
        data_item.smooth(window)
        data.append(data_item)
        
        cores.append(core)
        sections.append(section)
        faces.append(face)
        ACorDCs.append(ACorDC)


#%% define plotting function
def plotquarter(y_vec,ycor,d,meas,button,axs,rescale):
    
    width = y_vec[1] - y_vec[0]
    
    for y in y_vec:
        
        
        idx = ycor==y
        
        tmeas = meas[idx]
        tbut = button[idx]
        tycor = ycor[idx]
        td = d[idx]
        
        for i in range(len(tmeas)-1):
            
            if tbut[i] == 0 or sum(button)>0.9*len(button):
                
                axs.add_patch(Rectangle((y-(width-0.2)/2,td[i]),(width-0.2),td[i+1]-td[i],facecolor=my_cmap(rescale(tmeas[i]))))
            else:
                axs.add_patch(Rectangle((y-(width-0.2)/2,td[i]),(width-0.2),td[i+1]-td[i],facecolor='k'))
            
    return()

#%% define function to find unique elements in list
def unique(list1):
 
    # initialize a null list
    unique_list = []
 
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    
    return(unique_list)

#%% Make colormap

# make colormap
my_cmap = matplotlib.colormaps['Spectral']




#%% plot each section



for data_face in data:

    print('ALHIC2302 - '+data_face.section+' - '+ data_face.ACorDC+
          ' - '+data_face.face+' - '+str(window)+' mm smooth')
    
                
    idx1 = data_face.button_s == 0
    idx2 = data_face.y_s > data_face.y_vec[0]
    idx3 = data_face.y_s < data_face.y_vec[-1]
    idx = (idx1 * idx2) * idx3
    
    if sum(idx)>0:
        pltmin = np.percentile(data_face.meas_s[idx],5)
        pltmax = np.percentile(data_face.meas_s[idx],95)
    else:
        pltmin = np.percentile(data_face.meas_s,5)
        pltmax = np.percentile(data_face.meas_s,95)
    rescale = lambda k: (k-pltmin) /  (pltmax-pltmin)


   
    # make figure
    fig, ax = plt.subplots(figsize=(5,8),dpi=200)
    
    ax.set_xlim([min(data_face.y_vec)-5,max(data_face.y_vec)+5])
    ax.set_ylim([min(data_face.depth),max(data_face.depth)])


        
    # plot data
    plotquarter(data_face.y_vec,
                data_face.y_s,
                data_face.depth_s,
                data_face.meas_s,
                data_face.button_s,
                ax,
                rescale)
    
    # housekeeping
    ax.set_title('ALHIC2302 - '+data_face.section+' - '+data_face.ACorDC+' - '+data_face.face+' - '+str(window)+' mm smooth')




    # ad colorbar
    norm = matplotlib.colors.Normalize(vmin=pltmin,vmax=pltmax)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=my_cmap),ax=ax,
                  orientation='horizontal',label='Current (amps)')

    
    # save figure
    fname = path_to_figures +'alhic2302-'+data_face.section+'-'+data_face.ACorDC+'-'+data_face.face+'.png'
    fig.savefig(fname,bbox_inches='tight')
    plt.close(fig)

