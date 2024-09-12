#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 09:39:20 2024

Plot all data from ALHIC2302 Shallow Cores in an annimation for use in presentations

@author: Liam
"""

#%% Import packages 

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# general
import numpy as np
import pandas as pd
import os

# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# annimation
import moviepy.video.io.ImageSequenceClip

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

# progress bar
from tqdm import tqdm


#%% User Inputs

# smoothing window, mm
window = 10

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/annimation_2302/'
metadata_file = 'metadata.csv'

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
    if core == 'alhic2302' and row['idx_abs'] < 48 and row['idx_abs']>0 and row['ACorDC']=='AC':

        
        section = row['section']
        face = row['face']
        ACorDC = row['ACorDC']
        
        print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
    
        data_item = ECM(core,section,face,ACorDC)
        data_item.rem_ends(10)
        data_item.smooth(window)
        data_item.norm_outside()
        data.append(data_item)
        
        cores.append(core)
        sections.append(section)
        faces.append(face)
        ACorDCs.append(ACorDC)

sec = set(sections)

#%% define plotting function
def plotquarter(y_vec,ycor,d,meas,button,axs,rescale):
    
    width = y_vec[1] - y_vec[0]
    
    for y in y_vec:
        
        
        idx = ycor==y
        
        tmeas = meas[idx]
        tbut = button[idx]
        tycor = ycor[idx]
        td = d[idx]
        
        if False:
            for i in range(len(tmeas)-1):
                if tbut[i] == 0:
                    axs.add_patch(Rectangle((y-(width-0.2)/2,td[i]),(width-0.4),td[i+1]-td[i],facecolor=my_cmap(rescale(tmeas[i]))))
        else:
            for i in range(int(round(len(tmeas)/10)-1)):
                if tbut[i] == 0:
                    axs.add_patch(Rectangle((y-(width-0.2)/2,
                                             td[i*10]),(width-0.4),td[i*10+10]-td[i*10],
                                            facecolor=my_cmap(rescale(np.mean(tmeas[i*10:i*10+9])))))
            
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

#%% One big plot

if True:
    AC_all = []
    DC_all = []
    dmin = 100
    dmax = 0 
    for d in data:
        if d.core=='alhic2302':
            if d.ACorDC == 'AC':
                AC_all.extend(d.meas)
            else:
                DC_all.extend(d.meas)
            if min(d.depth)<dmin:
                dmin = min(d.depth)
            if max(d.depth)>dmax:
                dmax = max(d.depth)
    
    ACpltmin_all = np.percentile(AC_all,5)
    ACpltmax_all = np.percentile(AC_all,95)
    ACrescale_all = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)

    fig2,ax2 = plt.subplots(1, 2, 
                            gridspec_kw={'width_ratios': [3, 3]},
                            figsize=(6,6),
                            dpi=200)
    
    #fig2.suptitle('ALHIC2302 - All Data - '+str(window)+' mm smooth')
    ax2[0].set_title('AC - Left')
    ax2[1].set_title('AC - Top')

    
    # left-specific
    ax2[0].set_xlim([120, 0])
    ax2[1].yaxis.tick_right()
    ax2[1].yaxis.set_label_position("right")
    ax2[1].set_xlim([0,120])
    
    
    # applies to all
    for a in [ax2[0],ax2[1]]:
        a.set_ylabel('Depth (m)')
        a.set_xlabel('Distance From Center (mm)')




#%% plot each section

for sec in unique(sections):

    
    # print update
    print("Running Section "+sec)
    
    # set data to empty
    AC_t = None
    AC_l = None
    #loop through data 
    for d in data:
        
        # find faces
        if d.core=='alhic2302':
            if d.section==sec:
                if d.ACorDC == 'AC':
                    if d.face == 't':
                        AC_t = d
                    if d.face == 'l':
                        AC_l = d                    
    
    # find depth max and minimum
    minvec = []
    maxvec = []
    AC_all = []
    for data_face in [AC_t,AC_l]:
        if data_face != None:
            minvec.append(min(data_face.depth))
            maxvec.append(max(data_face.depth))
            
            AC_all.extend(data_face.meas)

    ACpltmin = np.percentile(AC_all,5)
    ACpltmax = np.percentile(AC_all,95)
    ACrescale = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)
    dmin = min(minvec)
    dmax = max(maxvec)
    
    
    if True:
        for a,data_face in zip([ax2[1],ax2[0]],[AC_t,AC_l]):
            
            if data_face != None:
                if data_face.face == 'l':
                    yall = data_face.y_right - data_face.y_s
                    yvec = data_face.y_right -  data_face.y_vec
                else:
                    yall = data_face.y_s -  data_face.y_left
                    yvec =data_face.y_vec -  data_face.y_left
                

                rescale = ACrescale_all

            
            
                # plot data
                plotquarter(yvec,
                            yall,
                            data_face.depth_s,
                            data_face.meas_s,
                            data_face.button_s,
                            a,
                            rescale)
        print("     done with big plot")
        
#%% Housekeeping 

fig2.suptitle('ALHIC2302')

fig2.tight_layout()
fig2.subplots_adjust(wspace=0)
ax2[0].set_ylim([20+2/2, 20-2/2])

#ACcbar_ax = fig2.add_axes([0.1,0.01,0.8,0.02])
fig2.subplots_adjust(bottom=0.22)
ACcbar_ax = fig2.add_axes([0.1,0.07,0.8,0.04])
ACnorm = matplotlib.colors.Normalize(vmin=ACpltmin_all,vmax=ACpltmax_all)
ACcbar = fig2.colorbar(matplotlib.cm.ScalarMappable(norm=ACnorm, cmap=my_cmap),cax=ACcbar_ax,
              orientation='horizontal',label='Current (amps)')



#%% Loop through

plt_depths = np.linspace(11.4,44,1600)
width = 2
dcnt=0
image_names=[]

pbar = tqdm(total=len(plt_depths), position=0, leave=True)


for d in plt_depths:

    
    
    ax2[0].set_ylim([d+width/2, d-width/2])
    ax2[1].set_ylim([d+width/2, d-width/2])
    ax2[0].yaxis.set_major_locator(plt.MultipleLocator(0.2))
    fig2.suptitle(f'ALHIC2302 AC ECM: {d:04.2f} m depth')
    

    fname = path_to_figures +'alhic2302_annimation_'+str(dcnt)+'.png'
    image_names.append(fname)
    fig2.savefig(fname)

    dcnt+=1
    pbar.update()
    
pbar.close()
    
#%% Make movie

image_files = [os.path.join(path_to_figures,img)
               for img in image_names]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files)
clip.write_videofile(path_to_figures+'/alhic2302_movie.mp4')

    