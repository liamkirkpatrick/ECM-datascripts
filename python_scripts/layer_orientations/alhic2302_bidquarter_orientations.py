#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 10:29:33 2024

New approach to calculating angle orientations on ALHIC quarter-core sections.

This approach first calcualtes a rough angle (within 1 degree) and 
then scans a narrower angle range with more accuracy


@author: Liam
"""

#%% Import packages

# general
import numpy as np
import pandas as pd
import math
from scipy.stats import pearsonr


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
window = 30

angle_res = 10

# spacing of line that we interpolate onto
interp_int = 0.0001

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
    face = row['face']
    
    #if core == 'alhic2302':
    if core == 'alhic2302' and (face == 'r' or face == 'tr'):

        
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



#%% Interpolate onto depth grid

def interp_meas(shifted_arr,depth,slopes,meas,y_vec,y_list):
    
    # make evenly spaced grid
    grid_min = math.ceil(min(depth)*100)/100
    grid_max = math.floor(max(depth)*100)/100
    grid = np.linspace(grid_min,
            grid_min + math.floor((grid_max-grid_min)/interp_int)*interp_int,
            math.floor((grid_max-grid_min)/interp_int)+1)
    
    # loop through angles
    meas_interp=[]
    scnt=0
    for s in slopes:
        
        # make empty numpy array
        temp = np.empty((len(grid),len(y_vec)))
        
        # loop through y_vec
        ycnt = 0
        for y in y_vec:
            idx = y_list==y
            
            temp[:,ycnt] = np.interp(grid,np.flip(shifted_arr[scnt][idx]),np.flip(meas[idx]))
            ycnt+=1
        
        meas_interp.append(temp)
        
        # count wwhich slope we're on
        scnt+=1
        
    return(meas_interp,grid)

#%% Get scores

def getscores(grid,meas_interp,slopes,dmin,dmax,y_vec):
    
    score_over_slopes = np.zeros([len(slopes),int((len(y_vec)-1)*len(y_vec)/2)])
    
    # loop through all slopes
    for s_cnt in range(len(slopes)):
        
        # make empty dataframe
        df_corr = pd.DataFrame()
        df_corr.drop(df_corr.index, inplace=True)
        
        # get appropriate array of measurments based on current angle
        array = meas_interp[s_cnt]
        
        # add those measurements to dataframe
        for i in range(len(y_vec)):
            idx1 = grid>=dmin
            idx2 = grid<=dmax
            df_corr[y_vec[i]] = array[idx1*idx2,i].tolist()
        
        # run multivariate correlations
        corr = df_corr.corr()
        corr = corr.to_numpy()
        n,m = np.shape(corr)
        
        # save all values from corner of corr array
        score_over_slopes[s_cnt,:] = corr[np.triu_indices(n,k=1)]
    
    return(score_over_slopes)


# #%% Calculate Angle

# def calc_angle(data_face,anglemin,anglemax,angle_res):
    
#     # make empty depth/angle 
#     score_list = []
#     angle=[]
#     angle_depth=[]
    
#     # assign angles to test
#     test_angle = np.linspace(anglemin,anglemax,int((anglemax-anglemin)*angle_res+1))
#     slopes = np.tan(test_angle * np.pi/180)
    
#     # assign local variables
#     depth = data_face.depth_s
#     meas = data_face.meas_s
#     y_list = data_face.y_s
#     y_vec = data_face.y_vec
    
#     # calculate shifts
#     shifts = np.tan(test_angle*np.pi/180)
    
#     for y in y_vec[0:]
    
#     # get full set of depths
#     # rounds min up to nearest whole cm, rounds max down to nearest whole cm
    
#     dmin = math.ceil((min(depth)+
#             (max(y_vec)-min(y_vec))/2/1000 * np.tan(max(test_angle)*np.pi/180)
#             )* 100)/100
#     dmax = math.ceil((max(depth)-
#             (max(y_vec)-min(y_vec))/2/1000 * np.tan(max(test_angle)*np.pi/180)
#             )* 100)/100
        
    
#     # make list of shifted arrays
#     shifted_arr = make_shifted_arrays(slopes,y_vec,depth,y_list)

#     # interpolate onto hifted arrays
#     meas_interp,grid = interp_meas(shifted_arr,depth,slopes,meas,y_vec,y_list)
    
#     # calculate scores
#     score_over_slopes = getscores(grid,meas_interp,slopes,dmin,dmax,y_vec)
    
#     # get min/max angle with max corr
#     max_indices = np.argmax(score_over_slopes, axis=0)
#     angles = test_angle[max_indices]
#     scores = np.max(score_over_slopes, axis=0)
    

    return angles,scores

# find longest stretch of zeros
def longest_stretch_zeros(data):
    max_length = 0
    current_length = 0
    start_index = 0
    best_start = 0
    best_end = 0

    for i, value in enumerate(data):
        if value == 0:
            if current_length == 0:
                start_index = i  # Start new sequence
            current_length += 1
        else:
            if current_length > max_length:
                max_length = current_length
                best_start = start_index
                best_end = start_index + current_length - 1
            current_length = 0

    # Check last sequence in case the longest sequence ends at the last element
    if current_length > max_length:
        best_start = start_index
        best_end = start_index + current_length - 1

    return list(range(best_start, best_end + 1))

#%% Compute Face angles

angles = []
scores = []

for data_face in data:
    
    print("Calculating Angle: "+data_face.core+
          ", section "+data_face.section+
          '-'+data_face.face+'-'+data_face.ACorDC)

    
    # assign angles to test
    anglemin = -75
    anglemax = 75
    test_angle = np.linspace(anglemin,anglemax,int((anglemax-anglemin)*angle_res+1))
    slopes = np.tan(test_angle * np.pi/180)
    
    # assign local variables
    depth = data_face.depth_s
    button = data_face.button_s
    meas = data_face.meas_s
    y_list = data_face.y_s
    y_vec = data_face.y_vec
    
    # make empty depth/angle 
    score_list = np.zeros(len(y_vec)-1)
    angle_list = np.zeros(len(y_vec)-1)
    
    # calculate shifts
    shifts = np.tan(test_angle*np.pi/180)
    
    # loop through each y track
    for i in range(len(y_vec[0:-1])):
        
        print("    Running tracks "+str(i)+" and "+str(i+1))
        
        # compute offset for each tracks
        offset =  shifts * ((y_vec[i+1]-y_vec[i])/2)/1000
        
        # find tracks
        track1_d = depth[y_list == y_vec[i]]
        track2_d = depth[y_list == y_vec[i+1]]
        track1_b = button[y_list == y_vec[i]]
        track2_b = button[y_list == y_vec[i+1]]
        track1_m = meas[y_list == y_vec[i]]
        track2_m= meas[y_list == y_vec[i+1]]
        
        # first, find longest stretch without button pressed.
        button_comb = track1_b+track2_b
        didx = longest_stretch_zeros(button_comb)
        
        
        # next, find dmin and dmax
        dmin = math.ceil( (min(track1_d[didx]) + max(offset)) * 100)/100
        dmax = math.ceil( (max(track1_d[didx]) - max(offset)) * 100)/100
        
        # only consider longer than 25cm
        if dmax-dmin > 0.15:
        
            # make grid to interpolate onto
            grid_min = math.ceil(min(depth)*100)/100
            grid_max = math.floor(max(depth)*100)/100
            grid = np.linspace(grid_min,
                    grid_min + math.floor((grid_max-grid_min)/interp_int)*interp_int,
                    math.floor((grid_max-grid_min)/interp_int)+1)
            
            for j in range(len(test_angle)):
            
                track1_doffset = track1_d + offset[j]
                track2_doffset = track2_d - offset[j]
            
                # interpolate onto grid
                track1_measinterp = np.interp(grid[(grid>=dmin)*(grid<=dmax)],
                                              np.flip(track1_doffset),
                                              np.flip(track1_m))
                track2_measinterp = np.interp(grid[(grid>=dmin)*(grid<=dmax)],
                                              np.flip(track2_doffset),
                                              np.flip(track2_m))
    
                # compute correlation
                rho, _ = pearsonr(track1_measinterp, track2_measinterp)   
    
                    
                if rho > score_list[i]:
                    score_list[i] = rho
                    angle_list[i] = test_angle[j]
        else:
            score_list[i] = np.nan
            angle_list[i] = np.nan
                
        print("    Angle is "+str(angle_list[i]))
    scores.append(score_list)
    angles.append(angle_list)

#%% compute angle for each face
face_angle = []
for a,s in zip(angles,scores):
    
    face_angle.append(np.nansum(a*s**2)/np.nansum(s**2))
    
#%% Test plots
fig, axs = plt.subplots(1,4,figsize=(12,9),dpi=100)

fig.suptitle('ALHIC2302 Layer Orientations')

# for all axis
for ax in axs:
    ax.grid()
    ax.set_ylabel('Depth (m)')
    ax.set_ylim([138,126])
    ax.set_xlim([-90,90])
axs[3].set_xlim([0,360])
axs[2].set_xlim(40,80)
    
axs[0].set_title('Right Face')
axs[1].set_title('Top Face')

cnt = 0
for data_face in data:
    
    if data_face.face == 'tr':
        ax = axs[1]
    else:
        ax = axs[0]
        
    if data_face.ACorDC == 'AC':
        c = 'r.'
        label = 'AC'
    else:
        c = 'b.'
        label = 'DC'
        
    mid_depth = (max(data_face.depth) + min(data_face.depth))/2
    angle = angles[cnt]
    score = scores[cnt]
    
    for i in range(len(angle)):
        

        ax.plot(angle[i],mid_depth,c,label=label,markersize=20*(score[i])**3)

    cnt+=1
        



handles, labels = axs[1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[1].legend(by_label.values(), by_label.keys())

handles, labels = axs[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[0].legend(by_label.values(), by_label.keys())

# axs[2].set_title('True Dip')
# axs[2].plot(true_AC_angle,true_AC_angledepth,'r.',label='AC')
# axs[2].plot(true_DC_angle,true_DC_angledepth,'b*',label='DC')
# axs[2].legend()

# axs[3].set_title('Orientation')
# axs[3].plot(true_AC_strike,true_AC_angledepth,'r.',label='AC')
# axs[3].plot(true_DC_strike,true_DC_angledepth,'b*',label='DC')
# axs[3].legend()

plt.tight_layout()


