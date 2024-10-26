"""


Angle Calculations on ALHIC2302

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

# weighted percentile statistics
from weighted_percentile import weighted_percentile

# math
from statsmodels.stats.weightstats import DescrStatsW

#%% User Inputs

# smoothing window, mm
window = 20

# paths
path_to_data = '../../data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/deep_orientations/'
metadata_file = 'metadata.csv'

# Correlation
# distance over which to complete correlation comparison
comp_range = 0.25

# depth interval to complete correlation comparison
comp_int = 0.001

score_floor = 0.65

# spacing of line that we interpolate onto
interp_int = 0.00025

# do I want plots
plot = True

# set sections to run
to_run = ['155_2','156_2','158','159_3']


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
    
    if core == 'alhic2302':

        
        section = row['section']
        face = row['face']
        ACorDC = row['ACorDC']
        
        data_item = ECM(core,section,face,ACorDC)
        
        if section in to_run and ACorDC == 'AC' and (face == 'r' or face == 'tr'):
            print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
            
            data_item.rem_ends(20)
            data_item.smooth(window)
            data.append(data_item)
            
            cores.append(core)
            sections.append(section)
            faces.append(face)
            ACorDCs.append(ACorDC)

sec = set(sections)

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

def compute_dip_angles(data,sections,core):
#if True:

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
        angle_res = 0.5
        angle_low = -60
        angle_high = 60
        test_angle = np.linspace(angle_low,angle_high,int((angle_high-angle_low)*angle_res+1))
        
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
        
            
        idx = np.where(~np.isnan(corr_coef[0,:]))[0]
        # find max angle for each pair of tracks
        icnt = 0
        angle=[]
        score = []
        for i in idx:
            maxidx = np.argmax(corr_coef[:,i])
            angle.append(test_angle[maxidx])
            score.append(corr_coef[maxidx,i])
    
            icnt+=1
        
        # add data to dataframe
        sec_idx = unique_sec.index(d.section)
        df.loc[sec_idx,'section'] = d.section
        df.loc[sec_idx,'depth'] = np.mean(d.depth)
        angle_rowname = d.ACorDC + '-'+d.face + '-angles'
        score_rowname = d.ACorDC + '-'+d.face + '-scores'
        length_rowname = d.ACorDC + '-'+d.face + '-length'
        if angle_rowname not in df.columns:
            df[angle_rowname] = None
            df[score_rowname] = None
            df[length_rowname] = None
        df.at[sec_idx, angle_rowname] = angle
        df.at[sec_idx, score_rowname] = score
        df.at[sec_idx, length_rowname] = length
        
        
        # calculate percentage spread between 10 and 90th percentile
        meas10 = np.percentile(d.meas[d.button==0],10)
        meas90 = np.percentile(d.meas[d.button==0],90)
        spread = meas90/meas10
        df.at[sec_idx,'10-90 percentile spread'] = spread
        
        # increment counter
        dcnt+=1
    

    return(df)

#%% Run 2302
core = 'alhic2302'
df = compute_dip_angles(data,sections,'alhic2302')

#%% Calculate True Angles

# make empty columns for true angles
for n in ['AC-true-angles','AC-true-scores',
   'AC-true-orientations','AC-mean-angle',
   'AC-tr-mean-angle','AC-r-mean-angle']:
    df[n]=None
    
side = 'r'

# loop through all sections and compute drue dip
for index, row in df.iterrows():
    
    #print("Running row "+str(row['section']))
    
    # only do AC (copied code so also built for DC - whoops!)
    for ACorDC in ['AC']:
        
        top_angle = row[ACorDC+'-tr-angles']
        top_score = row[ACorDC+'-tr-scores']
        side_angle = row[ACorDC+'-'+side+'-angles']
        side_score = row[ACorDC+'-'+side+'-scores']
        top_length = row[ACorDC+'-tr-length']
        side_length = row[ACorDC+'-'+side+'-length']
        
        #%% save mean angle from each face
        for angle,score,col in zip([top_angle,side_angle],[top_score,side_score],['AC-tr-mean-angle','AC-r-mean-angle']):
            weighted_stats = DescrStatsW(angle, weights=score, ddof=0)
            df.at[index,col] = weighted_stats.mean
        
        dip=[]
        orientation=[]
        score=[]
        for i in range(len(top_angle)):
            a1 = top_angle[i]
            for j in range(len(side_angle)):
                a2=side_angle[j]
                score.append((top_score[i]**2)*(side_score[j]**2)*top_length[i]*side_length[j])
                eps = np.arctan( 1/np.tan(a1*np.pi/180) * np.tan(a2 * np.pi/180)) * 180/np.pi
                true = np.arctan(np.tan(a1 * np.pi/180) / np.sin((90-eps)*np.pi/180))* 180/np.pi
                eps = eps+90
                if true<0:
                    eps = eps+180
                    true = true * -1
                orientation.append(eps)
                dip.append(true)
                
        # compute central estimate
        weighted_stats = DescrStatsW(dip, weights=score, ddof=0)
        df.at[index,ACorDC+'-mean-angle'] = weighted_stats.mean

        # save to dataframe
        df.at[index,ACorDC+'-true-angles'] = dip
        df.at[index,ACorDC+'-true-scores'] = score
        df.at[index,ACorDC+'-true-orientations'] = orientation
        
        

# save dataframe
df.to_pickle('../../data/angles/'+core+'deep_angles.df')

# save a small array with just the mean side and true angles
df2 = df[['section','depth','AC-mean-angle', 'AC-tr-mean-angle',
'AC-r-mean-angle']]
df2.to_csv('../../data/angles/'+core+'_deepangles_means.csv',index=False)


