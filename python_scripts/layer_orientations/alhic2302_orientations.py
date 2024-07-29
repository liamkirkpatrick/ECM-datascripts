"""


Angle Calculations on ALHIC2302

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

#%% User Inputs

# smoothing window, mm
window = 20

# paths
path_to_data = '../../data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/orientations/'
metadata_file = 'metadata.csv'

# Correlation
# distance over which to complete correlation comparison
comp_range = 0.25

# depth interval to complete correlation comparison
comp_int = 0.001

score_floor = 0.75

# spacing of line that we interpolate onto
interp_int = 0.00025

# do I want plots
plot = True

# set angles to cycle through
ang = 75
ang_res = 5 #fraction of angle
test_angle = np.linspace(-ang,ang,2*ang*ang_res+1)
slopes = np.tan(test_angle * np.pi/180)

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
        
        if max(data_item.depth) < 47 and min(data_item.depth)>3:
            print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
            
            data_item.rem_ends(10)
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

#%% set up for angle calcs

right_AC_angle = []
right_DC_angle = []
top_AC_angle = []
top_DC_angle = []

right_AC_angledepth = []
right_DC_angledepth = []
top_AC_angledepth = []
top_DC_angledepth = []

right_AC_anglescore = []
right_DC_anglescore = []
top_AC_anglescore = []
top_DC_anglescore = []



#%% loop through all data and calculate angles

# loop through each data file
for data_face in data:
    
    # print summary
    print("Calculating Angle: "+data_face.core+
          ", section "+data_face.section+
          '-'+data_face.face+'-'+data_face.ACorDC)
    
    # make empty depth/angle 
    score_list = []
    angle=[]
    angle_depth=[]
    
    # assign local variables
    depth = data_face.depth_s
    meas = data_face.meas_s
    y_list = data_face.y_s
    y_vec = data_face.y_vec
    
    # get full set of depths
    # rounds min up to nearest whole mm, rounds max down to nearest whole mm
    min_depth = math.ceil((min(depth) + 
                      (comp_range/2 + (max(y_vec)-min(y_vec))/2/1000
                       *np.tan(max(test_angle)*np.pi/180) + 0.0005))*100)/100
    max_depth = math.floor((max(depth) - 
                      (comp_range/2 + (max(y_vec)-min(y_vec))/2/1000
                       *np.tan(max(test_angle)*np.pi/180) + 0.0005))*100)/100
    print("min depth = "+str(min_depth)+", max depth = "+str(max_depth))
    print("Core min is: "+str(round(min(depth),3))+" and core max is: "+
          str(round(max(depth),3)))
    
    if min_depth>=max_depth:
        print("  *stick is too short*")
    else:
        
        # let's make our array of shifted values: 
        print("Making array of shifted Values")
        
        # make empty list of np arrays
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
            
        # now let's interpolate our shifted values onto evenly spaced grid
        print("Interpolating")
        
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
            
        # now let's loop through all of our depths and calculate a correlation score
        print("Calculating Score")
        
        # set depths to loop through
        d_vec = np.linspace(min_depth,
                        min_depth + math.floor((max_depth-min_depth)/comp_int)*comp_int,
                        math.floor((max_depth-min_depth)/comp_int)+1)
        d_vec = np.round(d_vec,6)
        
        # start progress bar
        pbar = tqdm(total=len(d_vec), position=0, leave=True)
        
        # Loop though depths
        for d in d_vec:

            dmax = d + comp_range/2
            dmin = d - comp_range/2
            
            score_over_slopes = []
            
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
                
                # count score (tril pulls out bottom left triangle)
                score_over_slopes.append(np.sum(np.sum(np.tril(corr))) / np.count_nonzero(np.tril(corr)))
                
            
            
            # add list of scores, angles, depths
            score_list.append(max(score_over_slopes))
            angle.append(test_angle[score_over_slopes.index(max(score_over_slopes))].tolist())
            angle_depth.append(d)
            
            # update progress bar
            pbar.update()
            
        complete = True
        pbar.close()
            
    # save in the right place
    if data_face.core=='alhic2302':
        if data_face.ACorDC == 'AC':
            if data_face.face == 't':
                top_AC_angle.extend(angle)
                top_AC_angledepth.extend(angle_depth)
                top_AC_anglescore.extend(score_list)
            if data_face.face == 'l':
                right_AC_angle.extend(angle)
                right_AC_angledepth.extend(angle_depth)
                right_AC_anglescore.extend(score_list)              
        else:
            if data_face.face == 't':
                top_DC_angle.extend(angle)
                top_DC_angledepth.extend(angle_depth)
                top_DC_anglescore.extend(score_list)
            if data_face.face == 'l':
                right_DC_angle.extend(angle)
                right_DC_angledepth.extend(angle_depth)
                right_DC_anglescore.extend(score_list)



    print("Complete")
    print("")


#%% Compute true angles

def find_index(float_list, target_float):
    try:
        index = float_list.index(target_float)
        return index
    except ValueError:
        return None

true_AC_angledepth = []
true_AC_angle = []
true_AC_strike = []

for d,a2,s in zip(top_AC_angledepth,top_AC_angle,top_AC_anglescore):
    
    if s>0.9:
    
        idx = find_index(right_AC_angledepth,d)
        
        if idx != None and right_AC_anglescore[idx]>score_floor and s>score_floor:
            
            # assign depth
            true_AC_angledepth.append(d)
            
            # get a2
            a1 = right_AC_angle[idx]
            
            # compute angle between a1 plane and true dip plane
            eps = np.arctan( 1/np.tan(a1*np.pi/180) * np.tan(a2 * np.pi/180)) * 180/np.pi
            
            # calc true dip
            true = np.arctan(np.tan(a1 * np.pi/180) / np.sin((90-eps)*np.pi/180))* 180/np.pi
           
            if true<0:
                eps = eps+180
                true = true * -1
            
            true_AC_angle.append(true)
            true_AC_strike.append(eps)
            

            
true_DC_angledepth = []
true_DC_angle = []
true_DC_strike = []

for d,a2,s in zip(top_DC_angledepth,top_DC_angle,top_DC_anglescore):
    
    if s>0.9:
    
        idx = find_index(right_DC_angledepth,d)
        
        if idx != None and right_DC_anglescore[idx]>score_floor and s>score_floor:
            
            # assign depth
            true_DC_angledepth.append(d)
            
            # get a2
            a1 = right_DC_angle[idx]
            
            # compute angle between a1 plane and true dip plane
            eps = np.arctan( 1/np.tan(a1*np.pi/180) * np.tan(a2 * np.pi/180)) * 180/np.pi
            
            # calc true dip
            true = np.arctan(np.tan(a1 * np.pi/180) / np.sin((90-eps)*np.pi/180))* 180/np.pi
            
            if true<0:
                eps = eps+180
                true = true * -1
            
            true_DC_angle.append(true)
            true_DC_strike.append(eps)
            

#%% Plot all

fig, axs = plt.subplots(1,4,figsize=(12,9),dpi=100)

fig.suptitle('ALHIC2302 Layer Orientations')

# for all axis
for ax in axs:
    ax.grid()
    ax.set_ylabel('Depth (m)')
    ax.set_ylim([45,15])
    ax.set_xlim([min(test_angle),max(test_angle)])

axs[3].set_xlim([0,360])
axs[2].set_xlim(40,80)
    
axs[0].set_title('Left Face')
for i in range(len(right_AC_angle)):
    if right_AC_anglescore[i] >score_floor:
        axs[0].plot(right_AC_angle[i],right_AC_angledepth[i],'r.',label='AC - high confidence')
    else:
        axs[0].plot(right_AC_angle[i],right_AC_angledepth[i],'k*',label='AC - poor confidence')
for i in range(len(right_DC_angle)):
    if right_DC_anglescore[i] >score_floor:
        axs[0].plot(right_DC_angle[i],right_DC_angledepth[i],'b.',label='DC - high confidence')
    else:
        axs[0].plot(right_DC_angle[i],right_DC_angledepth[i],'k*',label='DC - poor confidence')

axs[1].set_title('Top Face')
for i in range(len(top_AC_angle)):
    if top_AC_anglescore[i] >score_floor:
        axs[1].plot(top_AC_angle[i],top_AC_angledepth[i],'r.',label='AC - high confidence')
    else:
        axs[1].plot(top_AC_angle[i],top_AC_angledepth[i],'k*',label='poor confidence')
for i in range(len(top_DC_angle)):
    if top_DC_anglescore[i] >score_floor:
        axs[1].plot(top_DC_angle[i],top_DC_angledepth[i],'b.',label='DC - high confidence')
    else:
        axs[1].plot(top_DC_angle[i],top_DC_angledepth[i],'k*',label='DC - poor confidence')
        



handles, labels = axs[1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[1].legend(by_label.values(), by_label.keys())

handles, labels = axs[0].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[0].legend(by_label.values(), by_label.keys())

axs[2].set_title('True Dip')
axs[2].plot(true_AC_angle,true_AC_angledepth,'r.',label='AC')
axs[2].plot(true_DC_angle,true_DC_angledepth,'b*',label='DC')
axs[2].legend()

axs[3].set_title('Orientation')
axs[3].plot(true_AC_strike,true_AC_angledepth,'r.',label='AC')
axs[3].plot(true_DC_strike,true_DC_angledepth,'b*',label='DC')
axs[3].legend()



plt.tight_layout()

fig.savefig(path_to_figures+'2023_firstangleplot.png')



