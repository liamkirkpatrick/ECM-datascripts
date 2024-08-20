#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 19:27:25 2024

Goal here is to create a plot which walks readers through the layer 
orientation calculation proccess. Some experimentation expected as
to which subplots make the most sense and help build to the goal.

Will build mostly on master_orientations.py. This code is super awkward as a
result, but it makes sense to just hack bits off of that rather than building
something from the ground up.

@author: Liam
"""

#%% Import packages

# general
import numpy as np
import pandas as pd
import math
from scipy import stats

# progress bar
from tqdm import tqdm

# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

#%% user Inputs

# paths
path_to_data = '../../data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/methods_paper/'
metadata_file = 'metadata.csv'

# smoothing window
window = 10

#%% set up colormap

# make colormap
cmap = matplotlib.colormaps.get_cmap('coolwarm')

#%% Read in metadata and import data - ALHIC2302

meta = pd.read_csv(path_to_data+metadata_file)

# import each script as an ECM class item
data = []
cores = []
sections = []
faces = []
ACorDCs = []
max_tracks = 0
for index,row in meta.iterrows():
    
    core = row['core']
    section = row['section']
    section_num = section.split("_")
    face = row['face']
    face = row['face']
    ACorDC = row['ACorDC']
    
    #if core == 'alhic2302' and section=='51_2' and row['ACorDC']=='AC' and face=='l':
    if core == 'alhic2302' and int(section_num[0])==18 and ACorDC == 'AC' and face == 't':
        
        
        
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
        
        if len(data_item.y_vec)>max_tracks:
            max_tracks = len(data_item.y_vec)

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

#%% Make shifted array

def shift_depths(slopes,depth,y_vec,y_list):
    
    shifted_depth = []
    
    for s in slopes:
        slope_shifted = []
        
        
        for y in y_vec:
            idx = y_list==y
            slope_shifted.append(depth[idx]-(y-(y_vec[-1]+y_vec[0])/2)*s/1000)
        
        shifted_depth.append(slope_shifted)
        
    return(shifted_depth)

#%% Find longest stretch of 0s

def find_longest_zeros(button):
    max_length = 0
    current_length = 0
    start_index = 0

    best_start = 0
    best_end = 0

    for i, value in enumerate(button):
        if value == 0:
            if current_length == 0:
                start_index = i
            current_length += 1
        else:
            if current_length > max_length:
                max_length = current_length
                best_start = start_index
                best_end = i
            current_length = 0

    # Check the last segment in case the longest stretch of 0s ends at the last element
    if current_length > max_length:
        max_length = current_length
        best_start = start_index
        best_end = len(button)

    return list(range(best_start, best_end))

#%% Interpolate onto shifted arrays

def interp_onto_shifted(shifted_depth,depth,slopes,meas,y_vec,y_list,interp_int,button):
    
    # make evenly spaced grid
    grid_min = math.ceil(min(depth)*100)/100
    grid_max = math.floor(max(depth)*100)/100
    grid = np.linspace(grid_min,
            grid_min + math.floor((grid_max-grid_min)/interp_int)*interp_int,
            math.floor((grid_max-grid_min)/interp_int)+1)
    
    
    interp_meas = []
    interp_depth = []
    depth_min = np.ones(len(y_vec)) * 0
    depth_max = np.ones(len(y_vec)) * 2000
    
    
    # loop through all tracks
    scnt = 0
    for s in slopes:
                
        slope_shifted = shifted_depth[scnt]
        
        
        # loop through all tracks
        interp_meas_slope = []
        interp_depth_slope = []
        ycnt=0
        for y in y_vec:
            
            # find track index
            idx = y_list==y
            
            # pull out shifted depths from this track
            track_shifted = slope_shifted[ycnt]
            
            # filter for depths between button push
            but_idx = find_longest_zeros(button[idx])

            
            if len(but_idx)>0:

                # make grid
                grid_min = math.ceil(min(track_shifted[but_idx])*100)/100
                grid_max = math.floor(max(track_shifted[but_idx])*100)/100
                grid = np.linspace(grid_min,
                                   grid_min + math.floor((grid_max-grid_min)/interp_int)*interp_int,
                                   math.floor((grid_max-grid_min)/interp_int)+1)
                
                interp_meas_slope.append(np.interp(grid,np.flip(track_shifted),np.flip(meas[idx])))
                interp_depth_slope.append(grid)
                
                if grid_min >= depth_min[ycnt]:
                    depth_min[ycnt] = grid_min
                if grid_max <= depth_max[ycnt]:
                    depth_max[ycnt] = grid_max
            else:
                interp_meas_slope.append([np.NaN])
                interp_depth_slope.append([np.NaN])
            
            ycnt+=1
        
        interp_depth.append(interp_depth_slope)
        interp_meas.append(interp_meas_slope)
        scnt+=1
            
    return(interp_meas,interp_depth,depth_min,depth_max)

#%% Calculate angles


#if True:

column_names = ['section']
unique_sec = unique(sections)
df = pd.DataFrame(unique_sec,columns=column_names)
df['depth'] = 0

# loop through all data
dcnt=0
d = data[0]
    
# calculate angle
print("Calculating Angle - "+d.core+
      " - "+d.section+
      " - "+d.face+
      " - "+d.ACorDC)

length = []
angle_res = 10
angle_low = -75
angle_high = 75
test_angle = np.linspace(angle_low,angle_high,int(angle_high-angle_low*angle_res+1))

# assign angles to test
slopes = np.tan(test_angle * np.pi/180)

# assign local variables
depth = d.depth_s
meas = d.meas_s
button = d.button_s
y_list = d.y_s
y_vec = d.y_vec

# get shifted depths for each track
print("    getting shifted depth")
shifted_depth = shift_depths(slopes,depth,y_vec,y_list)

# get interpolated arrays
print("    getting interpolated arrays")
interp_int = 0.00025
interp_meas,interp_depth,depth_min,depth_max = interp_onto_shifted(shifted_depth,depth,
                                                slopes,meas,y_vec,y_list,
                                                interp_int,button)

# calc dip angle
print("    calc dip angle")
corr_coef=np.zeros((len(slopes),len(y_vec)**2)) * np.NaN
track_combos = []
scnt = 0
for s in slopes:
    
    interp_meas_slope = interp_meas[scnt]
    interp_depth_slope = interp_depth[scnt]
    
    
    for i in range(len(y_vec)):
        for j in range(len(y_vec)):
            
            # check we're not comparing a track with itself and 
            # check each track has reasonable length (>15cm) (being over 100m indicates it's a track of all button)
            if (i!=j 
                and (depth_max[i]-depth_min[i])>0.15
                and (depth_max[i]-depth_min[i])<2
                and (depth_max[j]-depth_min[j])>0.15
                and (depth_max[j]-depth_min[j])<2):
                
                if scnt==0:
                    track_combos.append('Tracks '+str(i)+' and '+str(j))
                
                length.append(depth_max[i]-depth_min[i])
                
                dmin = max([depth_min[i],depth_min[j]])+0.0010001
                dmax = min([depth_max[i],depth_max[j]])-0.0009999
                
                k = 1

                track1_idx = (interp_depth_slope[i]>=dmin)* (interp_depth_slope[i]<=dmax)
                track2_idx = (interp_depth_slope[j]>=dmin)* (interp_depth_slope[j]<=dmax)
                track1_meas = interp_meas_slope[i]
                track2_meas = interp_meas_slope[j]
                
                # track1_depth_all = interp_depth_slope[i]
                # track2_depth_all = interp_depth_slope[j]
                # track1_depth = track1_depth_all[track1_idx]
                # track2_depth = track2_depth_all[track2_idx]
                
                if np.isnan(np.sum(track1_meas[track1_idx])) or np.isnan(np.sum(track2_meas[track2_idx])):
                    print('       ************* Tracks'+str(i)+' and '+str(j)+ ' includes nans **************')
                elif sum(track1_idx)<10 or sum(track2_idx)<10:
                    print('       ************* No Significant Overlap ********************  ')
                else:
                    result = stats.pearsonr(track1_meas[track1_idx],track2_meas[track2_idx])
                    corr_coef[scnt,len(y_vec)*i+j] = result.statistic
    scnt+=1

    
idx_corr = np.where(~np.isnan(corr_coef[0,:]))[0]
# find max angle for each pair of tracks
icnt = 0
angle=[]
score = []
for i in idx_corr:
    maxidx = np.argmax(corr_coef[:,i])
    angle.append(test_angle[maxidx])
    score.append(corr_coef[maxidx,i])

    icnt+=1
    
    
#%% calc angle
angle = np.array(angle)
score = np.array(score)
true_angle = np.sum(angle*score)/np.sum(score)
pair_angle_idx = np.argmax(corr_coef[:,5])
pair_angle = test_angle[pair_angle_idx]
        

#%% Make Plot

fig,axs = plt.subplots(2,2,figsize=(10,10),dpi=200)
#fig,axs = plt.subplots(1,2,figsize=(10,25))
fig.suptitle('Demo Plot')


# suplot 1 - make plot
ycnt = 0
for y in d.y_vec:
    idx = d.y_s == y
    axs[0,0].plot(d.meas_s[idx],d.depth_s[idx],label="track "+str(ycnt+1),color=cmap(ycnt/len(y_vec)))
    ycnt+=1
b_t = 16.4
b_b = 16.65
b_l = 1.9e-8
b_r = 2.8e-8
axs[0,0].plot([b_l,b_l,b_r,b_r,b_l],[b_b,b_t,b_t,b_b,b_b],'k')
axs[0,0].set_ylim([np.max(d.depth),np.min(d.depth)])

# subplot 2 - make plot
axs[0,1].set_xlim([b_l,b_r])
axs[0,1].set_ylim([b_b,b_t])
idx = d.y_s == d.y_vec[0]
axs[0,1].plot(d.meas_s[idx],d.depth_s[idx],label="Track 1 ",color=cmap(0.0))
idx = d.y_s == d.y_vec[-1]
axs[0,1].plot(d.meas_s[idx],d.depth_s[idx],label="Track "+str(len(d.y_vec)),color=cmap(1.0))
shift_1 = np.tan(pair_angle*np.pi/180) * (d.y_left-d.y_vec[0])/1000
idx = d.y_s == d.y_vec[0]
axs[0,1].plot(d.meas_s[idx],d.depth_s[idx]+shift_1,label="Track 1 ",color=cmap(0.0),linestyle='dashed')
shift_6 = np.tan(pair_angle*np.pi/180) * (d.y_left-d.y_vec[-1])/1000
idx = d.y_s == d.y_vec[-1]
axs[0,1].plot(d.meas_s[idx],d.depth_s[idx]+shift_6,label="Track "+str(len(d.y_vec))+" shifted",color=cmap(1.0),linestyle='dashed')
    
# subplot 3 - make plot
for i in idx_corr:
    axs[1,0].plot(test_angle,corr_coef[:,i],'k-')
    axs[1,0].fill_between(test_angle,corr_coef[:,i],corr_coef[:,i]-2,color='0.1',alpha=0.03)
axs[1,0].plot(test_angle,corr_coef[:,5],'r-')

# figure settup
#       subplot 1
axs[0,0].set_xlabel('Conductivity (amps)')
axs[0,0].set_ylabel('depth')
axs[0,0].legend()
#       subplot 2
axs[0,1].set_xlabel('Conductivity')
axs[0,1].set_ylabel('depth')
axs[0,1].legend()
#       subplot 3
axs[1,0].set_xlabel('Test Angle (degrees)')
axs[1,0].set_ylabel('Correlation Coefficent')
axs[1,0].set_ylim([0,1.1])

plt.tight_layout()


# line between subplots
transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(axs[0,0].transData.transform([b_r,b_t]))
coord2 = transFigure.transform(axs[0,1].transData.transform([b_l,b_t,]))
coord3 = transFigure.transform(axs[0,0].transData.transform([b_r,b_b]))
coord4 = transFigure.transform(axs[0,1].transData.transform([b_l,b_b,]))
line1 = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                                transform=fig.transFigure,color='k')
line2 = matplotlib.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),
                                transform=fig.transFigure,color='k')

fig.lines.append(line1)
fig.lines.append(line2)






