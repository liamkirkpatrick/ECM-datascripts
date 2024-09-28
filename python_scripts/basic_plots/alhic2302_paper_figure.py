#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:39:20 2024

Modification of ALHIC2302 plotting to make a figure for 
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
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/three_dim/'
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
    
    section = row['section']
    section_num = section.split("_")
    
    # filter for ALHIC2302 shallow ice
    if core == 'alhic2302' and int(section_num[0]) < 30 and int(section_num[0])>26:

        
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
                axs.add_patch(Rectangle((y-(width-0.2)/2,td[i]),(width-0.2),td[i+1]-td[i],facecolor='w'))
            
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

if onebigplot:
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
    DCpltmin_all = np.percentile(DC_all,5)
    DCpltmax_all = np.percentile(DC_all,95)
    ACrescale_all = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)
    DCrescale_all = lambda k: (k-DCpltmin) /  (DCpltmax-DCpltmin)

    fig2,ax2 = plt.subplots(1, 5, 
                            gridspec_kw={'width_ratios': [3, 3,2, 3, 3]},
                            figsize=(8,15),
                            dpi=240)
    
    #fig2.suptitle('ALHIC2302 - All Data - '+str(window)+' mm smooth')
    ax2[2].axis('off')
    ax2[0].set_title('AC - Top')
    ax2[1].set_title('AC - Left')
    ax2[3].set_title('DC - Top')
    ax2[4].set_title('DC - Left')
    
    # left-specific
    for a in [ax2[0],ax2[3]]:
        #a.yaxis.tick_left()
        a.set_xlim([120, 0])
        
    # top specific
    for a in [ax2[1],ax2[4]]:
        a.yaxis.tick_right()
        a.yaxis.set_label_position("right")
        a.set_xlim([0,120])
        
    # applies to all
    for a in [ax2[0],ax2[1],ax2[3],ax2[4]]:
        if a is not ax2[1]:
            a.set_ylabel('Depth (m)')
        a.set_xlabel('Distance From Center (mm)')
        a.set_ylim([dmax, dmin])
        a.yaxis.set_major_locator(plt.MultipleLocator(0.2))



#%% plot each section



for sec in unique(sections):

    
    # print update
    print("Running Section "+sec)
    
    # set data to empty
    AC_t = None
    AC_l = None
    DC_t = None
    DC_l = None
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
                else:
                    if d.face == 't':
                        DC_t = d
                    if d.face == 'l':
                        DC_l = d
    
    # find depth max and minimum
    minvec = []
    maxvec = []
    AC_all = []
    DC_all = []
    for data_face in [AC_t,AC_l,DC_t,DC_l]:
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
        fig, ax = plt.subplots(1, 5, gridspec_kw={'width_ratios': [3, 3,2, 3, 3]},figsize=(9,6),dpi=200)
        
        # top-specific
        for a in [ax[0],ax[3]]:
            #a.yaxis.tick_right()
            a.set_xlim([120, 0])
            
        # right specific
        for a in [ax[1],ax[4]]:
            a.yaxis.tick_right()
            a.yaxis.set_label_position("right")
            a.set_xlim([0,120])
            
        # applies to all
        for a in [ax[0],ax[1],ax[3],ax[4]]:
            a.set_ylabel('Depth (m)')
            a.set_xlabel('Distance From Center (mm)',fontsize=6)
            a.set_ylim([dmax, dmin])
            
            
        for a,data_face in zip([ax[1],ax[0],ax[4],ax[3]],[AC_l,AC_t,DC_l,DC_t]):
            
            if data_face != None:
                if data_face.face == 'l':
                    yall = data_face.y_right - data_face.y_s
                    yvec = data_face.y_right -  data_face.y_vec
                else:
                    yall = data_face.y_s -  data_face.y_left
                    yvec =data_face.y_vec -  data_face.y_left
                
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
        ax[1].set_title('AC - Top')
        ax[0].set_title('AC - Left')
        ax[4].set_title('DC - Top')
        ax[3].set_title('DC - Left')
        
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
                      orientation='horizontal',label='Current (amps)')
        DCcbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=DCnorm, cmap=my_cmap),cax=DCcbar_ax,
                      orientation='horizontal',label='Current (amps)')
        
        # save figure
        fname = path_to_figures +'alhic2302-'+sec+'.png'
        fig.savefig(fname,bbox_inches='tight')
        plt.close(fig)
        print("     done with small plot")
    
    if onebigplot:
        
        #[ax[1],ax[0],ax[4],ax[3]],[AC_l,AC_t,DC_l,DC_t]
        
        for a,data_face in zip([ax2[1],ax2[0],ax2[4],ax2[3]],[AC_l,AC_t,DC_l,DC_t]):
            
            
            if data_face != None:
                if data_face.face == 'l':
                    yall = data_face.y_right - data_face.y_s
                    yvec = data_face.y_right -  data_face.y_vec
                else:
                    yall = data_face.y_s -  data_face.y_left
                    yvec =data_face.y_vec -  data_face.y_left
                
                if data_face.ACorDC =='AC':
                    rescale = ACrescale_all
                else:
                    rescale = DCrescale_all
            
            
                # plot data
                plotquarter(yvec,
                            yall,
                            data_face.depth_s,
                            data_face.meas_s,
                            data_face.button_s,
                            a,
                            rescale)
        print("     done with big plot")

#%% save all plot

if onebigplot:
    
    fig2.tight_layout()
    fig2.subplots_adjust(wspace=0)
    
    ACcbar_ax = fig2.add_axes([0.08,-0.01,0.35,0.005])
    ACnorm = matplotlib.colors.Normalize(vmin=ACpltmin_all,vmax=ACpltmax_all)
    DCcbar_ax = fig2.add_axes([0.57,-0.01,0.35,0.005])
    DCnorm = matplotlib.colors.Normalize(vmin=DCpltmin_all,vmax=DCpltmax_all)
    ACcbar = fig2.colorbar(matplotlib.cm.ScalarMappable(norm=ACnorm, cmap=my_cmap),cax=ACcbar_ax,
                  orientation='horizontal',label='Current (amps)')
    DCcbar = fig2.colorbar(matplotlib.cm.ScalarMappable(norm=DCnorm, cmap=my_cmap),cax=DCcbar_ax,
                  orientation='horizontal',label='Current (amps)')
        
    print("saving big plot")
    fname = path_to_figures +'alhic2302_paperplot.png'
    fig2.savefig(fname)
    print("done saving big plot")
    
    