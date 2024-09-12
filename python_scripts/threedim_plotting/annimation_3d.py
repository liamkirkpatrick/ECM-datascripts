#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 08:18:10 2024

Make 3D ECM annimation

@author: Liam
"""

#%% Import packages 

# general
import numpy as np
import pandas as pd
import scipy as sp
import re

# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as ticker

import moviepy.video.io.ImageSequenceClip



# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM
import os

#%% User Inputs

# smoothing window, mm
window = 10

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/annimation_3d/'
metadata_file = 'metadata.csv'

# ONE BIG PLOT
onebigplot = True
indplots = False

#%% Set up colormaps

cmap = matplotlib.colormaps['Spectral']

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

#%% Read in metadata and import data

meta = pd.read_csv(path_to_data+metadata_file)

# get dimension bounds
min_depth = 1000
max_depth = 0
top_width = 0
l_width = 0

# import each script as an ECM class item
data = []
cores = []
sections = []
faces = []
ACorDCs = []
for index,row in meta.iterrows():
    
    core = row['core']
    ACorDC = row['ACorDC']
    section = row['section']
    section_num = section.split("_")
    
    # filter for ALHIC2302 shallow ice
    #if core == 'alhic2302' and int(section_num[0]) < 30 and int(section_num[0])>27 and ACorDC == 'AC':
    if core == 'alhic2302' and ACorDC == 'AC' and row['idx_abs'] < 48 and row['idx_abs']>0:
    #if core == 'alhic2302' and ACorDC == 'AC' and row['idx_abs'] < 30 and row['idx_abs']>20:
        
        section = row['section']
        face = row['face']
        
        
        print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
    
        data_item = ECM(core,section,face,ACorDC)
        data_item.rem_ends(10)
        data_item.smooth(window)
        data.append(data_item)
        
        cores.append(core)
        sections.append(section)
        faces.append(face)
        ACorDCs.append(ACorDC)
        
        #check depths
        if min(data_item.depth)<min_depth:
            min_depth= min(data_item.depth)
        if max(data_item.depth)>max_depth:
            max_depth= max(data_item.depth)
            
            
        if data_item.face == 't' and data_item.y_right-data_item.y_left > top_width:
            top_width = data_item.y_right-data_item.y_left
        if data_item.face == 'l' and data_item.y_right-data_item.y_left > l_width:
            l_width = data_item.y_right-data_item.y_left

sec = set(sections)

# convert dimensions to mm
min_depth = min_depth *1000
max_depth = max_depth * 1000

#%% Make plot - outline and admin




# make Plot - data

# get overall colorscale
AC_all = []
for d in data:
    AC_all.extend(d.meas)
ACpltmin = np.percentile(AC_all,10)
ACpltmax = np.percentile(AC_all,90)
rescale = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)


#%% loop through depths and make figures
hi = 41*1000
lo = 25*1000
depths = np.linspace(lo,hi,int((hi-lo)/10+1))


cnt = 0
for depth in depths:
    
    print("Making plot "+str(cnt)+" of "+str(len(depths)))
    
    min_depth = depth - 1000
    max_depth = depth + 1000
    
    
    
    
    fig = plt.figure(figsize=(30,30),num=1,dpi=120)
    axs = fig.add_subplot(1,1,1,projection='3d')
    
    # loop through all sections
    for d in data:
        
        if (max(d.depth)>min_depth/1000 and max(d.depth)<max_depth/1000) or (min(d.depth)<max_depth/1000 and min(d.depth)>min_depth/1000):
        
            print("    Running Section: "+d.section+" - face: "+d.face)
            
            ycor = d.y_s
            y_vec = d.y_vec
            meas = d.meas_s
            button = d.button_s
            depth = d.depth_s
            yspace = y_vec[1]-y_vec[0]
            
            # assign cross
            if d.face =='t':
                z = [0,0,0,0]
            else: 
                x = [0,0,0,0]
                
            # loop trhough tracks
            for j in range(len(y_vec)):
                
                # pull out values
                ytrack = y_vec[j]
                idx = ycor==ytrack
                tmeas = meas[idx]
                tbut = button[idx]
                tycor = ycor[idx]
                td = depth[idx]
                
                # interpolate onto much rougher depth scale to speed plotting
                int_lo = round(min(td),2)
                int_hi = round(max(td),2)
                depth_interp = np.linspace(int_lo,int_hi,int((int_hi-int_lo)/0.005)+1)
                meas_interp = np.interp(depth_interp,np.flip(td),np.flip(tmeas))
                
                # loop through all datapoints
                for i in range(len(meas_interp)-1):
                    
                    if depth_interp[i]*1000>min_depth and depth_interp[i]*1000<max_depth:
                    
                        # assign depth (y)values
                        y = [depth_interp[i]* 1000, depth_interp[i+1]* 1000, depth_interp[i+1]* 1000, depth_interp[i]* 1000] 
                        
                        # assign values in face direction
                        if d.face == 't':
                            
                            x = [(ytrack-d.y_left-yspace*0.45),
                                 (ytrack-d.y_left-yspace*0.45),
                                 (ytrack-d.y_left+yspace*0.45),
                                 (ytrack-d.y_left+yspace*0.45)]
                            
                        else:
                            
                            z = [-(d.y_right-ytrack-yspace*0.45),
                                 -(d.y_right-ytrack-yspace*0.45),
                                 -(d.y_right-ytrack+yspace*0.45),
                                 -(d.y_right-ytrack+yspace*0.45),
                                 ]
            
                        # make the plot
                        verts = [list(zip(x,y,z))]
                        collection = Poly3DCollection(verts, linewidths=1, alpha=1,antialiased=False)
                        face_color = cmap(rescale(meas_interp[i]))
                        collection.set_facecolor(face_color)
                        axs.add_collection3d(collection)
    

    # Make plot - outline
    
    # Make outline - length plot
    axs.plot([0,0],[min_depth,max_depth],[0,0],color='k')
    axs.plot([top_width,top_width],[min_depth,max_depth],[0,0],color='k')
    axs.plot([0,0],[min_depth,max_depth],[-l_width,-l_width],color='k')
    
    # make outline - crosspieces
    axs.plot([0,top_width],[min_depth,min_depth],[0,0],color='k')
    axs.plot([0,top_width],[max_depth,max_depth],[0,0],color='k')
    axs.plot([0,0],[min_depth,min_depth],[0,-l_width],color='k')
    axs.plot([0,0],[max_depth,max_depth],[0,-l_width],color='k')
    
    # set parameters 
    vert_dist = 50
    hor_dist = 10
    cent_offset_vert = 10
    cent_offset_hor = 20
    
    # calculate curve
    theta0 = 3*np.pi/2 - np.arctan(cent_offset_hor/(l_width-cent_offset_vert))
    radius = np.sqrt(cent_offset_hor**2 + (l_width-cent_offset_vert)**2)
    vert_dist = np.sqrt(radius**2 - (top_width-cent_offset_hor)**2) + cent_offset_vert
    #vert_dist = (top_width-cent_offset_hor) * np.tan(np.arccos((top_width-cent_offset_hor)/radius))+ cent_offset_vert
    theta1 = 2*np.pi - np.arccos((top_width-cent_offset_hor)/radius)
    #theta1 = 2*np.pi
    theta = np.linspace(theta0,theta1,40)
    
    # vertical outside line
    axs.plot([top_width,top_width],[min_depth,min_depth],[0,-vert_dist],color='k')
    axs.plot([top_width,top_width],[max_depth,max_depth],[0,-vert_dist],color='k')
    
    # curve
    x = radius * np.cos(theta)
    z = radius * np.sin(theta)
    axs.plot(x+cent_offset_hor, min_depth * np.ones_like(x), z-cent_offset_vert, color='k')
    axs.plot(x+cent_offset_hor, max_depth * np.ones_like(x), z-cent_offset_vert, color='k')
    
    # Make plot - final admin
    
    # turn off grid
    axs.grid(False)
    
    
    fsz = 30
    
    # x-axis labels
    hi = int(np.round(max_depth/1000)*1000)
    lo = int(np.round(min_depth/1000)*1000)
    ticks = np.linspace(lo,hi,int(np.round((hi-lo)/1000))+1)
    axs.set_yticks(ticks)
    #axs.set_yticks(np.round(ticks/1000,1),fontsi)
    axs.set_yticklabels(np.round(ticks/1000,1),fontsize=fsz)
    axs.tick_params(axis='y', pad=120)
    axs.set_ylabel('Depth (m)',fontsize=fsz,labelpad=380)
    
    axs.set_xticks([])
    axs.set_zticks([])
    
    axs.set_ylim([min_depth,max_depth])
    
    axs.view_init(elev=25, azim=120)

    # set to scale
    axs.set_aspect('equal')
    
    figfile = path_to_figures+'ann_3d_'+str(cnt)+'.png'
    fig.savefig(figfile, 
            transparent = False,  
            facecolor = 'white'
            )
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()
    
    cnt+=1
    
#%% make movie

def get_png_files(folder_path):
    return [file for file in os.listdir(folder_path) if file.endswith('.png')]

image_names = get_png_files(path_to_figures)
image_names = sorted(image_names, key=lambda x: int(re.findall(r'\d+', x)[-1]))
image_files = [os.path.join(path_to_figures,img)
               for img in image_names]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files,fps=24)
clip.write_videofile('threedim_movie.mp4',fps=24)

