#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:02:45 2024

Script to consolodate all sampling data into a readable format

@author: Liam
"""

#%% Import Packages

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


#%% User Inputs

# toggle messages on/off
debug = True

# filepaths
path_to_data = '../../data/sampling/'
path_to_angle = '../../data/angles/'

# Set datasets to include in master
discrete_datasets = ['water_iso']

#%% read in angles

core = 'alhic2302'
angles = pd.read_csv(path_to_angle+core+'_deepangles_means.csv',index_col=0)


#%% Do the Thing

# loop through all of the datasets
d = discrete_datasets[0]

# read in values
values = pd.read_csv(path_to_data+d+'/'+d+'_values.csv',index_col=0)
depths = pd.read_excel(path_to_data+d+'/'+d+'_depths.xlsx',sheet_name='depths',index_col=0)
meta = pd.read_excel(path_to_data+d+'/'+d+'_depths.xlsx',sheet_name='metadata',index_col=0)

# pick out sticks
#sticks = meta['stick'].to_list()

# loop through sticks
for s in meta.index:

    if debug:
        print("Running Stick: " + s)
    
    # pull out stick info
    parts = s.split('-')
    core = parts[0]
    section = parts[1]
    cut = parts[2]

    if debug:
        print("    Core: "+core)
        print("    Section: "+section)
        print("    Cut: "+cut)

    # pull out x and y dimension

    # assign angles
    if section in angles.index:

        run_angle = True
        x_angle = angles.at[section,'AC-tr-mean-angle']
        y_angle = angles.at[section,'AC-r-mean-angle']

        # get cut data
        y = (meta.at[s,'y_hi'] + meta.at[s,'y_lo'])/2
        x = (meta.at[s,'x_hi'] + meta.at[s,'x_lo'])/2

        if debug:
            print("    y_angle = "+str(y_angle))
            print("    x_angle = "+str(x_angle))
            print("    y = "+str(y))
            print("    x = "+str(x))

        # compute shift due to layer angle adjustment
        y_shift = -y * np.tan(y_angle*np.pi/180)
        x_shift = x * np.tan(x_angle*np.pi/180)
        shift = y_shift + x_shift
        if debug:
            print("    Y_shift = "+str(round(y_shift,3)))
            print("    X_shift = "+str(round(x_shift,3)))

    else:


        if section == '230_4':
            run_angle = True
            x_angle = -30
            y_angle = -30
            
            # get cut data
            y = (meta.at[s,'y_hi'] + meta.at[s,'y_lo'])/2
            x = (meta.at[s,'x_hi'] + meta.at[s,'x_lo'])/2

            # compute shift due to layer angle adjustment
            y_shift = -y * np.tan(y_angle*np.pi/180)
            x_shift = x * np.tan(x_angle*np.pi/180)
            shift = y_shift + x_shift
            if debug:
                print("    Y_shift = "+str(round(y_shift,3)))
                print("    X_shift = "+str(round(x_shift,3)))

        else:
            print("    ERROR")
            print("    No angle data for this section")
            run_angle = False

    
    td = depths[s+'_td']
    bd = depths[s+'_bd']
    td = td[td.notna()].to_list()
    bd = bd[bd.notna()].to_list()

    for i in range(len(td)):

        # add section info to sheet
        values.loc[s+'-'+str(i+1),'core'] = core
        values.loc[s+'-'+str(i+1),'section'] = section
        values.loc[s+'-'+str(i+1),'cut'] = cut
        values.loc[s+'-'+str(i+1),'sample_number'] = i+1
        values.loc[s+'-'+str(i+1),'sampleID'] = s+'-'+str(i+1)

        # add depths (not dip adjusted)
        values.loc[s+'-'+str(i+1),'top_depth'] = td[i]
        values.loc[s+'-'+str(i+1),'bottom_depth'] = bd[i]
        values.loc[s+'-'+str(i+1),'ave_depth'] = (td[i]+bd[i])/2

        # add depths (dip adjusted) and x/y dimension
        #print("    Run Angle = "+str(run_angle))
        if run_angle:
            values.loc[s+'-'+str(i+1),'top_depth_adj'] = td[i]+shift/1000
            values.loc[s+'-'+str(i+1),'bottom_depth_adj'] = bd[i]+shift/1000
            values.loc[s+'-'+str(i+1),'ave_depth_adj'] = (td[i]+bd[i])/2+shift/1000

            # save x and y dimension (unit meters)
            values.loc[s+'-'+str(i+1),'x_m'] = x/1000
            values.loc[s+'-'+str(i+1),'y_m'] = y/1000
        else:
            values.loc[s+'-'+str(i+1),'top_depth_adj'] = np.nan
            values.loc[s+'-'+str(i+1),'bottom_depth_adj'] = np.nan
            values.loc[s+'-'+str(i+1),'ave_depth_adj'] = np.nan
            values.loc[s+'-'+str(i+1),'x_m'] = np.nan
            values.loc[s+'-'+str(i+1),'y_m'] = np.nan


#%% 
# define non-linear deuterium excess
values['dln'] = 1000*(np.log(values['dD']/1000+1) - (8.47* np.log(values['d18O']/1000+1) - 28.5*(np.log(values['d18O']/1000+1))**2))

#%% Save output

# sort columns into desired order
order = ['core','section','cut','sample_number','sampleID','top_depth','bottom_depth','ave_depth','top_depth_adj','bottom_depth_adj','ave_depth_adj','x_m','y_m','dD','d18O','dxs','dln']
values = values[order]

# sort rows into desired order
values = values.sort_values(by=['core','section','cut','sample_number'])

# save
values.to_csv(path_to_data+'/'+d+'/master_'+d+'.csv')