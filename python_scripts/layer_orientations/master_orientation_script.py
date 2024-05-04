#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:05:44 2024

A more effective overall approach to orientations 

Orient#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/orientations/'
metadata_file = 'metadata.csv'

# smoothing window
window = 10



# spacing of line that we interpolate onto
interp_int = 0.00025

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
    face = row['face']
    
    #if core == 'alhic2302':
    if core == 'alhic2302' and int(section_num[0])<24 and int(section_num[0])>22:
        
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


#%% Calculate dip from interpolated array


def calc_dip(interp_depth,interp_meas,shifted_depth,slopes,y_vec,depth_min,depth_max):
    
    corr_coef=np.zeros((len(slopes),len(y_vec)**2)) * np.NaN
    
        
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
        
                    dmin = max([depth_min[i],depth_min[j]])+0.01001
                    dmax = min([depth_max[i],depth_max[j]])-0.00999
                    
                    k = 1

                    track1_idx = (interp_depth_slope[i]>=dmin)* (interp_depth_slope[i]<=dmax)
                    track2_idx = (interp_depth_slope[j]>=dmin)* (interp_depth_slope[j]<=dmax)
                    track1_meas = interp_meas_slope[i]
                    track2_meas = interp_meas_slope[j]
                    
                    # track1_depth_all = interp_depth_slope[i]
                    # track2_depth_all = interp_depth_slope[j]
                    # track1_depth = track1_depth_all[track1_idx]
                    # track2_depth = track2_depth_all[track2_idx]
                    
                    result = stats.pearsonr(track1_meas[track1_idx],track2_meas[track2_idx])
                    corr_coef[scnt,6*i+j] = result.statistic
        scnt+=1
    

        
    return(corr_coef)

    

#%% Calculate rough angle

def get_roughangle(data_face):
    

    
    # assign angles to test
    ang = 75
    test_angle = np.linspace(-ang,ang,2*ang+1)
    slopes = np.tan(test_angle * np.pi/180)
    
    # assign local variables
    depth = data_face.depth_s
    meas = data_face.meas_s
    button = data_face.button_s
    y_list = data_face.y_s
    y_vec = data_face.y_vec
    
    # get shifted depths for each track
    print("        getting shifted depth")
    shifted_depth = shift_depths(slopes,depth,y_vec,y_list)
    
    # get interpolated arrays
    print("        getting interpolated arrays")
    interp_meas,interp_depth,depth_min,depth_max = interp_onto_shifted(shifted_depth,depth,
                                                   slopes,meas,y_vec,y_list,
                                                   interp_int,button)
    
    # calc dip angle
    print("        calc dip angle")
    corr_coef = calc_dip(interp_depth,interp_meas,shifted_depth,slopes,y_vec,depth_min,depth_max)
    
    idx = np.where(~np.isnan(corr_coef[0,:]))[0]
    
    
    angle = []
    score = []
    for i in idx:
        
        maxidx = np.argmax(corr_coef[:,i])
        
        angle.append(test_angle[maxidx])
        score.append(corr_coef[maxidx,i])
        
    p20 = np.percentile(angle, 20)
    p80 = np.percentile(angle, 80)
    roughangle = np.sum(np.array(angle) * np.array(score)) / sum(score)
    
    

    return roughangle,p80,p20

#%% get real angle

def get_realangle(data_face,test_angle):
    

    
    # assign angles to test
    slopes = np.tan(test_angle * np.pi/180)
    
    # assign local variables
    depth = data_face.depth_s
    meas = data_face.meas_s
    button = data_face.button_s
    y_list = data_face.y_s
    y_vec = data_face.y_vec
    
    # get shifted depths for each track
    print("        getting shifted depth")
    shifted_depth = shift_depths(slopes,depth,y_vec,y_list)
    
    # get interpolated arrays
    print("        getting interpolated arrays")
    interp_meas,interp_depth,depth_min,depth_max = interp_onto_shifted(shifted_depth,depth,
                                                   slopes,meas,y_vec,y_list,
                                                   interp_int,button)
    
    # calc dip angle
    print("        calc dip angle")
    corr_coef = calc_dip(interp_depth,interp_meas,shifted_depth,slopes,y_vec,depth_min,depth_max)
    
    
    idx = np.where(~np.isnan(corr_coef[0,:]))[0]
    
    
    angle = []
    score = []
    for i in idx:
        
        maxidx = np.argmax(corr_coef[:,i])
        
        angle.append(test_angle[maxidx])
        score.append(corr_coef[maxidx,i])

    return angle,score
        
        
#%% core function
def calc_angles(data,sections,face,ACorDC):
    
    #df = pd.DataFrame(0)
    
    angles = []
    anglescores = []
    
    # loop through all data
    for d in data:
        
        # calculate angle
        print("Calculating Angle - "+d.core+
              " - "+d.section+
              " - "+d.face+
              " - "+d.ACorDC)
        
        # calculate rough angle
        print("    Getting Rough Angle")
        roughangle,anglemin,anglemax = get_roughangle(d)
        print("            angle = "+str(round(roughangle,1)))
        print("            angle min = "+str(anglemin))
        print("            angle max = "+str(anglemax))
        
        print("    Getting True Angle")
        angle_res = 10
        test_angle = np.linspace(round(anglemin),round(anglemax),int((round(anglemax)-round(anglemin))/angle_res+1))
        angles,scores = get_realangle(d,test_angle)
        
        
    df = 0
    return(df)

#%% test

df = calc_angles(data,sections,face,ACorDC)




