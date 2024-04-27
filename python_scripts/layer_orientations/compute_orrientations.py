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

angle_res = 10

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
    if core == 'alhic2302' and int(section_num[0])<54 and int(section_num[0])>20:
        
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

#%% Make shifted arrays

def make_shifted_arrays(slopes,y_vec,depth,y_list):
    
    shifted_arr = []
    
    # let's calcuate the shfited arrays. Loop through all slopes
    cnt = 0
    for s in slopes:
        
        temp = np.empty
        
        #loop through all tracks
        for y in y_vec:
            
            # find index of all points in this track
            idx = y_list==y
            
            # find all points for this track and add to the temp array,
            # after shifting by the appropriate depth for this slope
            if y==y_vec[0]:
                temp = depth[idx] - (y-(y_vec[-1]+y_vec[0])/2) * s/1000
            else:
                temp = np.append(temp,depth[idx] - (y-(y_vec[-1]+y_vec[0])/2) * s/1000)
            
        # now add this altered depth vector to our list of shifted arrays. 
        # This is in a debugging check for now
        if len(temp) == len(depth):
            shifted_arr.append(temp)
        else:
            print("Depths Don't Match!")
        
        # increment counter
        cnt+=1
        
    return shifted_arr

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





#%% Make function that determines the rough angle

def get_roughangle(data_face):
    
    # make empty depth/angle 
    score_list = []
    angle=[]
    angle_depth=[]
    
    # assign angles to test
    ang = 75
    test_angle = np.linspace(-ang,ang,2*ang+1)
    slopes = np.tan(test_angle * np.pi/180)
    
    # assign local variables
    depth = data_face.depth_s
    meas = data_face.meas_s
    y_list = data_face.y_s
    y_vec = data_face.y_vec
    
    # get full set of depths
    # rounds min up to nearest whole cm, rounds max down to nearest whole cm
    
    dmin = math.ceil((min(depth)+
            (max(y_vec)-min(y_vec))/2/1000 * np.tan(max(test_angle)*np.pi/180)
            )* 100)/100
    dmax = math.ceil((max(depth)-
            (max(y_vec)-min(y_vec))/2/1000 * np.tan(max(test_angle)*np.pi/180)
            )* 100)/100
        
    
    # make list of shifted arrays
    shifted_arr = make_shifted_arrays(slopes,y_vec,depth,y_list)

    # interpolate onto hifted arrays
    meas_interp,grid = interp_meas(shifted_arr,depth,slopes,meas,y_vec,y_list)
    
    # calculate scores
    score_over_slopes = getscores(grid,meas_interp,slopes,dmin,dmax,y_vec)
    
    # get min/max angle with max corr
    max_indices = np.argmax(score_over_slopes, axis=0)
    score_list = test_angle[max_indices]
    
    anglemin = min(score_list)
    anglemax = max(score_list)

    return anglemin,anglemax

#%% Calculate Angle

def calc_angle(data_face,anglemin,anglemax,angle_res):
    
    # make empty depth/angle 
    score_list = []
    angle=[]
    angle_depth=[]
    
    # assign angles to test
    ang = 75
    test_angle = np.linspace(anglemin,anglemax,int((anglemax-anglemin)*angle_res+1))
    slopes = np.tan(test_angle * np.pi/180)
    
    # assign local variables
    depth = data_face.depth_s
    meas = data_face.meas_s
    y_list = data_face.y_s
    y_vec = data_face.y_vec
    
    # get full set of depths
    # rounds min up to nearest whole cm, rounds max down to nearest whole cm
    
    dmin = math.ceil((min(depth)+
            (max(y_vec)-min(y_vec))/2/1000 * np.tan(max(test_angle)*np.pi/180)
            )* 100)/100
    dmax = math.ceil((max(depth)-
            (max(y_vec)-min(y_vec))/2/1000 * np.tan(max(test_angle)*np.pi/180)
            )* 100)/100
        
    
    # make list of shifted arrays
    shifted_arr = make_shifted_arrays(slopes,y_vec,depth,y_list)

    # interpolate onto hifted arrays
    meas_interp,grid = interp_meas(shifted_arr,depth,slopes,meas,y_vec,y_list)
    
    # calculate scores
    score_over_slopes = getscores(grid,meas_interp,slopes,dmin,dmax,y_vec)
    
    # get min/max angle with max corr
    max_indices = np.argmax(score_over_slopes, axis=0)
    angles = test_angle[max_indices]
    scores = np.max(score_over_slopes, axis=0)
    

    return angles,scores

#%% Compute top and side angles

angles = []
scores = []
faces = []
sections = []
meastypes = []
mid_depths = []

for data_face in data:
    
    print("Calculating Angle: "+data_face.core+
          ", section "+data_face.section+
          '-'+data_face.face+'-'+data_face.ACorDC)

    anglemin,anglemax = get_roughangle(data_face)


    angle,score= calc_angle(data_face,anglemin,anglemax,angle_res)
    
    angles.append(angle)
    scores.append(score)
    
    faces.append(data_face.face)
    sections.append(data_face.section)
    meastypes.append(data_face.ACorDC)
    
    mid_depths.append((max(data_face.depth) + min(data_face.depth))/2)

#%% Calculate True Dip

def calc_trueangle(top_angle,left_angle,top_score,left_score):
    
    angles = []
    scores = []
    epss = []
    
    for j in range(len(top_angle)):
        for k in range(len(left_angle)):
            
            eps = np.arctan( 1/np.tan(top_angle[j] *np.pi/180) * np.tan(left_angle[k] * np.pi/180)) * 180/np.pi
            angle = np.arctan(np.tan(top_angle[j] * np.pi/180) / np.sin((90-eps)*np.pi/180))* 180/np.pi
            
            score = top_score[j] * left_score[k]
            
            if angle<0:
                angle = angle * -1
                eps+=180
                
            if eps<0:
                eps+=360
            if eps>360:
                eps-=360
            
            epss.append(eps)
            angles.append(angle)
            scores.append(score)
    
    return epss,angles,scores


    
    return 
    
unique_section = unique(sections)

ac_true_angle = []
ac_true_score = []
dc_true_angle = []
dc_true_score = []
ac_true_eps = []
dc_true_eps = []
mid_depths = []

def find_matching_index(sections, faces, meastypes, sec2, face2, meas2):
    for index, (section, face, meastype) in enumerate(zip(sections, faces, meastypes)):
        if section == sec2 and face == face2 and meastype == meas2:
            return index
    return -1  # Return -1 if no match is found

for i in range(len(unique_section)):
    
    print("True Dip on section "+unique_section[i])
    
    sec_idx = (sections==unique_section[i])
    
    idx = find_matching_index(sections, faces, meastypes, unique_section[i], 't', 'AC')
    data_face = data[idx]
    mid_depths.append( (max(data_face.depth) + min(data_face.depth))/2)
    top_AC_score = scores[idx]
    top_AC_angle = angles[idx]
    
    idx = find_matching_index(sections, faces, meastypes, unique_section[i], 't', 'DC')
    top_DC_score = scores[idx]
    top_DC_angle = angles[idx]
    
    idx = find_matching_index(sections, faces, meastypes, unique_section[i], 'l', 'AC')
    left_AC_score = scores[idx]
    left_AC_angle = angles[idx]
    
    idx = find_matching_index(sections, faces, meastypes, unique_section[i], 'l', 'DC')
    left_DC_score = scores[idx]
    left_DC_angle = angles[idx]
    
    ac_eps,ac_angle,ac_score = calc_trueangle(top_AC_angle,left_AC_angle,top_AC_score,left_AC_score)
    ac_true_eps.append(ac_eps)
    ac_true_angle.append(ac_angle)
    ac_true_score.append(ac_score)
    
    dc_eps,dc_angle,dc_score = calc_trueangle(top_DC_angle,left_DC_angle,top_DC_score,left_DC_score)
    dc_true_eps.append(dc_eps)
    dc_true_angle.append(dc_angle)
    dc_true_score.append(dc_score)

#%% Test Plot

standard = 0.97
standard_true = standard**2

fig, axs = plt.subplots(1,4,figsize=(12,8),dpi=100)

fig.suptitle('ALHIC2302 Layer Orientations')

# for all axis
for ax in axs:
    ax.grid()
    ax.set_ylabel('Depth (m)')
    ax.set_ylim([45,15])
    ax.set_xlim([-90,90])
    
axs[0].set_title('Left Face')
axs[1].set_title('Top Face')
axs[0].set_xlabel('Aparent Dip Angle (degrees)')
axs[1].set_xlabel('Aparent Angle (degrees)')


cnt = 0
for data_face in data:
    
    if data_face.face == 't':
        ax = axs[1]
    else:
        ax = axs[0]
        
    if data_face.ACorDC == 'AC':
        c = 'r.'
        label = 'High Confidence (AC)'
        markersize=16
    else:
        c = 'b.'
        label = 'High Confidence (DC)'
        markersize = 10
        
    mid_depth = (max(data_face.depth) + min(data_face.depth))/2
    angle = angles[cnt]
    score = scores[cnt]
    
    for i in range(len(angle)):
        
        if score[i] >=standard:
            ax.plot(angle[i],mid_depth,c,label=label,markersize=markersize)
        else:
            ax.plot(angle[i],mid_depth,'k.',label='Low Confidence')
    cnt+=1
        

handles, labels = axs[1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[1].legend(by_label.values(), by_label.keys(),loc='upper center')

handles, labels = axs[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[0].legend(by_label.values(), by_label.keys(),loc='upper center')



axs[2].set_xlim([0,90])
axs[2].set_title('True Dip Angle')
axs[2].set_ylabel('Depth(m)')
axs[2].set_xlabel('True Dip Angle (degrees)')


for i in range(len(unique_section)):
    
    ac_angle = ac_true_angle[i]
    dc_angle = dc_true_angle[i]
    ac_score = ac_true_score[i]
    dc_score = dc_true_score[i]
    
    for j in range(len(ac_angle)):
        if ac_score[j] > standard_true:
            axs[2].plot(ac_angle[j],mid_depths[i],'r.',label='High Confidence (AC)',markersize=16)
        else:
            axs[2].plot(ac_angle[j],mid_depths[i],'k.',label='Low Confidence',markersize=2)
            
    for j in range(len(dc_angle)):
        if dc_score[j] > standard_true:
            axs[2].plot(dc_angle[j],mid_depths[i],'b.',label='High Confidence (DC)',markersize=10)
        else:
            axs[2].plot(dc_angle[j],mid_depths[i],'k.',label='Low Confidence',markersize=2)

handles, labels = axs[2].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[2].legend(by_label.values(), by_label.keys(),loc='upper center')


# Note where orientation is lost
o_lost = [19.73,20.39,24.28,30.72,39.53,40.47,41.39,42.7]

axs[3].set_title('Oreintation')
axs[3].set_xlim([-10,370])
axs[3].set_xlabel('Azimuth (degrees)')
axs[3].set_ylabel('Depth(m)')

for i in range(len(unique_section)):
    
    ac_eps = ac_true_eps[i]
    dc_eps = dc_true_eps[i]
    
    ac_angle = ac_true_angle[i]
    dc_angle = dc_true_angle[i]
    
    ac_score = ac_true_score[i]
    dc_score = dc_true_score[i]
    
    for j in range(len(ac_eps)):
        if ac_score[j] > standard_true:
            axs[3].plot(ac_eps[j],mid_depths[i],'r.',label='High Confidence (AC)',markersize=16)
        else:
            axs[3].plot(ac_eps[j],mid_depths[i],'k.',label='Low Confidence',markersize=2)
            
    for j in range(len(dc_eps)):
        if dc_score[j] > standard_true:
            axs[3].plot(dc_eps[j],mid_depths[i],'b.',label='High Confidence (DC)',markersize=10)
        else:
            axs[3].plot(dc_eps[j],mid_depths[i],'k.',label='Low Confidence',markersize=2)

for d in o_lost:
    axs[3].plot([-10, 370],[d,d],'g-',label='Orientation Lost (core log)',linewidth=3)


handles, labels = axs[3].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[3].legend(by_label.values(), by_label.keys(),loc='upper center')




plt.tight_layout()

fig.savefig(path_to_figures+'2023_firstangleplot.png')



