#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:39:20 2024

Plot all data from ALHCI2302 BID full rounds (quarter core measurements)

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
window = 20

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/first_plot_2302_BIDquarter/'
metadata_file = 'metadata.csv'

# ONE BIG PLOT
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
    section = row['section']
    face = row['face']
    ACorDC = row['ACorDC']
    
    if core == 'alhic2302' and (face == 'r' or face == 'tr') and 'repeat' not in section:
        
        print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
    
        data_item = ECM(core,section,face,ACorDC)
        data_item.rem_ends(10)
        data_item.smooth(window)
        data_item.norm_all()
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
        
        for i in range(len(tmeas)-1):
            
            if tbut[i] == 0:
                
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

for sec in unique(sections):

    
    # print update
    print("Running Section "+sec)
    
    # set data to empty
    AC_t = None
    AC_r = None
    DC_t = None
    DC_r = None
    #loop through data 
    for d in data:
        
        # find faces
        if d.core=='alhic2302':
            if d.section==sec:
                if d.ACorDC == 'AC':
                    if d.face == 'tr':
                        AC_t = d
                    if d.face == 'r':
                        AC_r = d                    
                else:
                    if d.face == 'tr':
                        DC_t = d
                    if d.face == 'r':
                        DC_r = d
    
    # find depth max and minimum
    minvec = []
    maxvec = []
    AC_all = []
    DC_all = []
    for data_face in [AC_t,AC_r,DC_t,DC_r]:
        if data_face != None:
            minvec.append(min(data_face.depth))
            maxvec.append(max(data_face.depth))
            
            if data_face.ACorDC == 'AC':
                AC_all.extend(data_face.meas)
            else:
                DC_all.extend(data_face.meas)
    ACpltmin = np.percentile(AC_all,5)
    ACpltmax = np.percentile(AC_all,95)
    DCpltmin = np.percentile(DC_all,5)
    DCpltmax = np.percentile(DC_all,95)  
    ACrescale = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)
    DCrescale = lambda k: (k-DCpltmin) /  (DCpltmax-DCpltmin)
    dmin = min(minvec)
    dmax = max(maxvec)
    
    if indplots:
   
        # make figure
        fig, ax = plt.subplots(1, 5, gridspec_kw={'width_ratios': [3, 3,3, 3, 3]},figsize=(9,6),dpi=200)
        
        # right-specific
        for a in [ax[0],ax[3]]:
            #a.yaxis.tick_right()
            a.set_xlim([120, 0])
            
        # top specific
        for a in [ax[1],ax[4]]:
            a.yaxis.tick_right()
            a.yaxis.set_label_position("right")
            a.set_xlim([0,120])
            
        # applies to all
        for a in [ax[0],ax[1],ax[3],ax[4]]:
            a.set_ylabel('Depth (m)')
            a.set_xlabel('Distance From Center (mm)',fontsize=6)
            a.set_ylim([dmax, dmin])
            
            
        for a,data_face in zip([ax[0],ax[1],ax[3],ax[4]],[AC_r,AC_t,DC_r,DC_t]):
            
            if data_face != None:
                if data_face.face == 'r':
                    yall = data_face.y_s - data_face.y_left
                    yvec = data_face.y_vec - data_face.y_left
                else:
                    yall = data_face.y_right - data_face.y_s
                    yvec = data_face.y_right - data_face.y_vec
                
                if data_face.ACorDC =='AC':
                    rescale = ACrescale
                else:
                    rescale = DCrescale
            
            
                # plot data
                plotquarter(yvec,
                            yall,
                            data_face.depth_s,
                            data_face.meas_s,
                            data_face.button_s,
                            a,
                            rescale)
        
        # housekeeping
        fig.suptitle('ALHIC2302 - '+sec+' - '+str(window)+' mm smooth')
        ax[2].axis('off')
        ax[0].set_title('AC - Right')
        ax[1].set_title('AC - Top')
        ax[3].set_title('DC - Right')
        ax[4].set_title('DC - Top')
        
        fig.tight_layout()
        plt.subplots_adjust(wspace=0)
    
        # ad colorbar
        #fig.subplots_adjust(bottom=0.8)
        #    ACcbar_ax = fig.add_axes([0.08,-0.07,0.35,0.05])
        ACcbar_ax = fig.add_axes([0.07,-0.05,0.35,0.05])
        ACnorm = matplotlib.colors.Normalize(vmin=ACpltmin,vmax=ACpltmax)
        DCcbar_ax = fig.add_axes([0.58,-0.05,0.35,0.05])
        DCnorm = matplotlib.colors.Normalize(vmin=DCpltmin,vmax=DCpltmax)
        ACcbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=ACnorm, cmap=my_cmap),cax=ACcbar_ax,
                      orientation='horizontal',label='Current (normalized)')
        DCcbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=DCnorm, cmap=my_cmap),cax=DCcbar_ax,
                      orientation='horizontal',label='Current (normalized)')
        
        # save figure
        fname = path_to_figures +'alhic2302-'+sec+'.png'
        #fig.savefig(fname,bbox_inches='tight')
        # plt.close(fig)
        
        plt.show()
    

    
    