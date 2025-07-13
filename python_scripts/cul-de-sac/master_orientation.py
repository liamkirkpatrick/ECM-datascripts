#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 11:23:00 2024

Orientation Script - try many
Doing things without functions - that seemed to be slowing things down

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
path_to_angles = '../../data/angles/'
metadata_file = 'metadata.csv'

# smoothing window
window = 10


#%% 
# define function to find unique elements in list
def unique(list1):
 
    # initialize a null list
    unique_list = []
 
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    
    return(unique_list)

#%% 
# Make shifted array

def shift_depths(slopes,depth,y_vec,y_list):
    
    shifted_depth = []
    
    for s in slopes:
        slope_shifted = []
        
        for y in y_vec:
            idx = y_list==y
            slope_shifted.append(depth[idx]-(y-(y_vec[-1]+y_vec[0])/2)*s/1000)
        
        shifted_depth.append(slope_shifted)
        
    return(shifted_depth)

#%% 
# Find longest stretch of 0s

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
                
                if grid_min >= depth_min[ycnt]: #YCNT ASSIGNMENT ERROR
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

def compute_dip_angles(data,sections,core):

    column_names = ['section']
    unique_sec = unique(sections)
    df = pd.DataFrame(unique_sec,columns=column_names)
    df['depth'] = 0
    
    # loop through all data
    dcnt=0
    for d in data:
        
        # calculate angle
        print("Calculating Angle - "+d.core+
              " - "+d.section+
              " - "+d.face+
              " - "+d.ACorDC)
        
        length = []
        angle_res = 0.1
        angle_low = -75
        angle_high = 75
        test_angle = np.linspace(angle_low,angle_high,int((angle_high-angle_low)/angle_res+1))
        
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
        interp_int = 0.001
        interp_meas,interp_depth,depth_min,depth_max = interp_onto_shifted(shifted_depth,depth,
                                                        slopes,meas,y_vec,y_list,
                                                        interp_int,button)
        
        # print depth_min and depth_max
        #print(depth_min)
        #print(depth_max)
        
        # calc dip angle
        print("    calc dip angle")
        corr_coef=np.zeros((len(slopes),len(y_vec)**2)) * np.NaN
        length_array=np.zeros((len(slopes),len(y_vec)**2)) * np.NaN
        y_offset_array=np.zeros((len(slopes),len(y_vec)**2)) * np.NaN
        scnt = 0
        for s in slopes:
            
            interp_meas_slope = interp_meas[scnt]
            interp_depth_slope = interp_depth[scnt]
            
            for i in range(len(y_vec)-1):
                for j in range(i+1,len(y_vec)):
                    
                    # check each track is not crazy long (all button is 100m)
                    if ((depth_max[i]-depth_min[i])<2
                        and (depth_max[j]-depth_min[j])<2 #):
                        and abs(y_vec[i]-y_vec[j]) > 5):
                        
                        if scnt==0 and False:
                            print('Tracks '+str(i)+' and '+str(j))
                        
                        #length.append(depth_max[i]-depth_min[i])
                        
                        dmin = max([depth_min[i],depth_min[j]])+0.0010001
                        dmax = min([depth_max[i],depth_max[j]])-0.0009999


                        # ensure there is significant overlap
                        if dmax-dmin>0.3:
                            
                            track1_idx = (interp_depth_slope[i]>=dmin)* (interp_depth_slope[i]<=dmax)
                            track2_idx = (interp_depth_slope[j]>=dmin)* (interp_depth_slope[j]<=dmax)
                            track1_meas = interp_meas_slope[i]
                            track2_meas = interp_meas_slope[j]
                        
                            
                            if np.isnan(np.sum(track1_meas[track1_idx])) or np.isnan(np.sum(track2_meas[track2_idx])):
                                print('       ************* Tracks'+str(i)+' and '+str(j)+ ' includes nans **************')
                            elif sum(track1_idx)<10 or sum(track2_idx)<10:
                                print('       ************* No Significant Overlap ********************  ')
                            else:
                                result = stats.pearsonr(track1_meas[track1_idx],track2_meas[track2_idx])
                                corr_coef[scnt,len(y_vec)*i+j] = result.statistic
                                length_array[scnt,len(y_vec)*i+j] = dmax - dmin #depth_max[i]-depth_min[i]
                                y_offset_array[scnt,len(y_vec)*i+j] = y_vec[j]-y_vec[i]
            scnt+=1
        
            
        idx = np.where(~np.isnan(corr_coef[0,:]))[0]
        # find max angle for each pair of tracks
        icnt = 0
        angle=[]
        score = []
        length = []
        y_offset = []

        for i in idx:
            maxidx = np.argmax(corr_coef[:,i])
            angle.append(test_angle[maxidx])
            score.append(corr_coef[maxidx,i])
            length.append(length_array[maxidx,i])
            y_offset.append(y_offset_array[maxidx,i])
    
            icnt+=1
        
        # add data to dataframe
        sec_idx = unique_sec.index(d.section)
        df.loc[sec_idx,'section'] = d.section
        df.loc[sec_idx,'depth'] = np.mean(d.depth)
        angle_rowname = d.ACorDC + '-'+d.face + '-angles'
        score_rowname = d.ACorDC + '-'+d.face + '-scores'
        length_rowname = d.ACorDC + '-'+d.face + '-length'
        y_offset_rowname = d.ACorDC + '-'+d.face + '-y_offset'
        if angle_rowname not in df.columns:
            df[angle_rowname] = None
            df[score_rowname] = None
            df[length_rowname] = None
            df[y_offset_rowname] = None
        df.at[sec_idx, angle_rowname] = angle
        df.at[sec_idx, score_rowname] = score
        df.at[sec_idx, length_rowname] = length
        df.at[sec_idx, y_offset_rowname] = y_offset
        
        # increment counter
        dcnt+=1
    
    df.to_pickle('../../data/angles/'+core+'_angles.df')

#%%
# Read in metadata and import data - ALHIC2201

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
    ACorDC = row['ACorDC']
    
    if core == 'alhic2401' and ACorDC == 'AC':
        
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
        
        if len(data_item.y_vec)>max_tracks:
            max_tracks = len(data_item.y_vec)

compute_dip_angles(data,sections,'alhic2401')

#%% 
# Read in metadata and import data - ALHIC2416

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
    
    if core == 'alhic2416':
        
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
        
        if len(data_item.y_vec)>max_tracks:
            max_tracks = len(data_item.y_vec)

compute_dip_angles(data,sections,'alhic2416')
