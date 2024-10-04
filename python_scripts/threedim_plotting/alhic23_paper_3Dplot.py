#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 08:18:10 2024

Make 3D ECM plot for paper

@author: Liam
"""

#%%
# Import packages 

# general
import numpy as np
import pandas as pd
import scipy as sp

# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as ticker

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
    if core == 'alhic2302' and int(section_num[0]) < 30 and int(section_num[0])>27 and ACorDC == 'AC':

        
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

fig = plt.figure(figsize=(30,30),num=1,dpi=200)
axs = fig.add_subplot(1,1,1,projection='3d')


# make Plot - data

# get overall colorscale
AC_all = []
for d in data:
    AC_all.extend(d.meas)
ACpltmin = np.percentile(AC_all,10)
ACpltmax = np.percentile(AC_all,90)
rescale = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)


if True:
    # loop through all sections
    for d in data:
        
        print("Running Section: "+d.section+" - face: "+d.face)
        
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
            depth_interp = np.linspace(int_lo,int_hi,int((int_hi-int_lo)/0.002)+1)
            meas_interp = np.interp(depth_interp,np.flip(td),np.flip(tmeas))
        
            
            # loop through all datapoints
            for i in range(len(meas_interp)-1):
                
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
            
            # # loop through all datapoints
            # for i in range(len(tmeas)-1):
                
            #     # assign depth (y)values
            #     y = [td[i]* 1000, td[i+1]* 1000, td[i+1]* 1000, td[i]* 1000] 
                
            #     # assign values in face direction
            #     if d.face == 't':
                    
            #         x = [(ytrack-d.y_left-yspace*0.45),
            #              (ytrack-d.y_left-yspace*0.45),
            #              (ytrack-d.y_left+yspace*0.45),
            #              (ytrack-d.y_left+yspace*0.45)]
                    
    
            #     else:
                    
            #         z = [-(d.y_right-ytrack-yspace*0.45),
            #              -(d.y_right-ytrack-yspace*0.45),
            #              -(d.y_right-ytrack+yspace*0.45),
            #              -(d.y_right-ytrack+yspace*0.45),
            #              ]
    
            #     # make the plot
            #     verts = [list(zip(x,y,z))]
            #     collection = Poly3DCollection(verts, linewidths=1, alpha=1,antialiased=False)
            #     face_color = cmap(rescale(tmeas[i]))
            #     collection.set_facecolor(face_color)
            #     axs.add_collection3d(collection)

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





fsz = 8

# x-axis labels
hi = int(np.round(max_depth/1000,1)*1000)
lo = int(np.round(min_depth/1000,1)*1000)
ticks = np.linspace(lo,hi,int(np.round((hi-lo)/200))+1)
axs.set_yticks(ticks)
#axs.set_yticks(np.round(ticks/1000,1),fontsi)
axs.set_yticklabels(np.round(ticks/1000,1),fontsize=fsz)
axs.tick_params(axis='y', pad=30)
axs.set_ylabel('Depth (m)',fontsize=fsz,labelpad=110)

axs.set_xticks([])
axs.set_zticks([])

axs.set_ylim([min_depth,max_depth])

# # z axis labels
# tick_spacing = 50
# axs.zaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
# ticks = np.round(axs.get_zticks(),2)
# axs.set_zticklabels(ticks,fontsize=fsz)
# axs.set_zlabel('  (mm)',fontsize=fsz,rotation=90,labelpad=15)
# axs.tick_params(axis='z', pad=10)

# # x axis labels
# tick_spacing = 50
# axs.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
# ticks = np.round(axs.get_xticks(),2)
# axs.set_xticklabels(ticks,fontsize=fsz,rotation=-15)
# axs.tick_params(axis='x', pad=5)
# axs.set_xlabel('(mm)',fontsize=fsz,rotation=90)#,labelpad=0)

#axs.set_xticklabels(ticks,fontsize=fsz,rotation=45)

axs.view_init(elev=25, azim=120)



# set to scale
axs.set_aspect('equal')

# set axis limits
#axs.set_zlim(-l_width,0)
#axs.set_xlim(0,top_width)

#fig.show()
figfile = path_to_figures+'alhic2302_3d.png'
fig.savefig(figfile, 
        transparent = False,  
        facecolor = 'white'
            )