#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plot all data from PICO2303 Firn Core

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

# import pywavelets
import pywt

#%% User Inputs

# smoothing window, mm
window = 10

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/first_plot_2302/'
metadata_file = 'metadata.csv'

# ONE BIG PLOT
onebigplot = True
indplots = False

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
    if core == 'pico2303':

        
        section = row['section']
        face = row['face']
        ACorDC = row['ACorDC']
        
        print("Reading "+core+", section "+str(section)+'-'+face+'-'+ACorDC)
    
        data_item = ECM(core,str(section),face,ACorDC)
        data_item.rem_ends(10)
        data_item.smooth(window)
        data.append(data_item)
        
        cores.append(core)
        sections.append(section)
        faces.append(face)
        ACorDCs.append(ACorDC)

sec = set(sections)

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


#%% Sort Data

DC_depth = []
DC_meas  = []

AC_depth_0 = []
AC_meas_0  = []
AC_depth_1 = []
AC_meas_1  = []
AC_depth_2 = []
AC_meas_2  = []


for d in data:
    
    
    if d.ACorDC == 'DC':
        
        for i in range(len(d.depth_s)):
            DC_depth.append(d.depth_s[i])
            DC_meas.append(d.meas_s[i])
            
        sorted_indices = sorted(range(len(DC_depth)), key=lambda k: DC_depth[k])
        DC_depth = [DC_depth[i] for i in sorted_indices]
        DC_meas = [DC_meas[i] for i in sorted_indices]
        
    else:
        
        idx1 = d.y_s == d.y_vec[1]
        mean1 = np.mean(d.meas_s[idx1])
        AC_depth_1.extend(list(d.depth_s[idx1]))
        AC_meas_1.extend(list(d.meas_s[idx1]))
        sorted_indices = sorted(range(len(AC_depth_1)), key=lambda k: AC_depth_1[k])
        AC_depth_1 = [AC_depth_1[i] for i in sorted_indices]
        AC_meas_1 = [AC_meas_1[i] for i in sorted_indices]
        
        idx0 = d.y_s == d.y_vec[0]
        mean0 = np.mean(d.meas_s[idx0])
        AC_depth_0.extend(list(d.depth_s[idx0]))
        AC_meas_0.extend(list(d.meas_s[idx0]/ mean0 * mean1))
        sorted_indices = sorted(range(len(AC_depth_0)), key=lambda k: AC_depth_0[k])
        AC_depth_0 = [AC_depth_0[i] for i in sorted_indices]
        AC_meas_0 = [AC_meas_0[i] for i in sorted_indices]

        idx2 = d.y_s == d.y_vec[2]
        mean2 = np.mean(d.meas_s[idx2])
        AC_depth_2.extend(list(d.depth_s[idx2]))
        AC_meas_2.extend(list(d.meas_s[idx2]/ mean2 * mean1))
        sorted_indices = sorted(range(len(AC_depth_2)), key=lambda k: AC_depth_2[k])
        AC_depth_2 = [AC_depth_2[i] for i in sorted_indices]
        AC_meas_2 = [AC_meas_2[i] for i in sorted_indices]
        
# sort





#%% Plot


fig,axs = plt.subplots(1,2,figsize=(10,25))

axs[0].set_xlabel('DC Conductivity')
axs[1].set_xlabel('AC Conductivity')

axs[0].plot(DC_meas,DC_depth,'k-')

axs[1].plot(AC_meas_0,AC_depth_0,'r-')
axs[1].plot(AC_meas_2,AC_depth_2,'b-')
axs[1].plot(AC_meas_1,AC_depth_1,'k-')


for ax in axs:
    ax.set_ylim([20.8,6.5])
    ax.set_ylim([20.8,6.5])
    ax.set_ylabel('Depth (m)')
fig.tight_layout()
fig.suptitle('PICO2303')
fig.subplots_adjust(top=0.97)
    

