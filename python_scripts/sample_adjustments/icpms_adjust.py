

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat October 9

@author: Liam
"""

#%% load packages

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%
# Set filepaths

path_to_data = '../../data/sampling/icpms/'
path_to_angle = '../../data/angles/'

cores = ['alhic2302','alhic2302','alhic2302','alhic2302','alhic1901']
sections = ['155_2','156_2','158','159_3','230_4']

#%%
# Load metadata
meta = pd.read_excel(path_to_data + 'icpms_metadata.xlsx')
meta['section'] = meta['section'].astype(str)

# load science data
df = pd.read_csv(path_to_data + 'icpms_raw.csv')
print(df.columns)
df['top_depth'] = df['Top Depth (cm)'] / 100
df['bot_depth'] = df['Bot Depth (cm)'] / 100
df['ave_depth'] = (df['Top Depth (cm)'] + df['Bot Depth (cm)']) / 200
df['y_m'] = 0
df['x_m'] = 0
df = df.drop(columns=['Top Depth (cm)','Bot Depth (cm)'])

# load angle data
core = 'alhic2302'
angles = pd.read_csv(path_to_angle+core+'_deepangles_means.csv',index_col=0)

# %%
# Apply dip adnjustment

# loop through sticks and core/section combos
for core,section in zip(cores,sections):

    # get the angle of the core, apply a dip adjustment if angle data exists
    if section in angles.index:

        run_angle = True
        x_angle = angles.at[section,'AC-tr-mean-angle']
        y_angle = angles.at[section,'AC-r-mean-angle']

        # get cut data
        meta_row = meta[(meta['core']==core) & (meta['section']==section)]
        y = (meta_row['y_hi'].to_numpy() + meta_row['y_lo'].to_numpy())/2
        x = (meta_row['x_hi'].to_numpy() + meta_row['x_lo'].to_numpy())/2
        top_depth = meta_row['top_depth'].to_numpy() + meta_row['offset_mm'].to_numpy()/1000
        #top_depth = meta_row['top_depth'].to_numpy() + meta_row['offset_mm'].to_numpy()/1000
        #bot_depth = meta_row['bot_depth'].to_numpy() + meta_row['offset_mm'].to_numpy()/1000

        #print(meta_row)

        # compute shift due to layer angle adjustment
        y_shift = -y * np.tan(y_angle*np.pi/180)
        x_shift = x * np.tan(x_angle*np.pi/180)
        shift = y_shift + x_shift

        print("Shift for section "+ section + " is "+str(shift))

    else:
        print("Angle not found for section "+ section)
        run_angle = False

        if section == '230_4':
            run_angle = True
            x_angle = -30
            y_angle = -30
            
            # get cut data
            meta_row = meta[(meta['core']==core) & (meta['section']==section) ]
            y = (meta_row['y_hi'].to_numpy() + meta_row['y_lo'].to_numpy())/2
            x = (meta_row['x_hi'].to_numpy() + meta_row['x_lo'].to_numpy())/2
            top_depth = meta_row['top_depth'].to_numpy() + meta_row['offset_mm'].to_numpy()/1000
            #top_depth = meta_row['top_depth'].to_numpy() + meta_row['offset_mm'].to_numpy()/1000
            #bot_depth = meta_row['bot_depth'].to_numpy() + meta_row['offset_mm'].to_numpy()/1000

            # compute shift due to layer angle adjustment
            y_shift = -y * np.tan(y_angle*np.pi/180)
            x_shift = x * np.tan(x_angle*np.pi/180)
            shift = y_shift + x_shift

            print("Shift for section "+ section + " is "+str(shift))

    df.loc[(df['core'] == core) & (df['section'] == section), 'ave_depth'] += top_depth
    if run_angle:
        df.loc[(df['core'] == core) & (df['section'] == section) , 'top_depth_adj'] += top_depth+shift/1000
        df.loc[(df['core'] == core) & (df['section'] == section), 'bot_depth_adj'] += top_depth+shift/1000
        df.loc[(df['core'] == core) & (df['section'] == section), 'ave_depth_adj'] += top_depth+shift/1000
        df.loc[(df['core'] == core) & (df['section'] == section), 'top_depth'] += top_depth
        df.loc[(df['core'] == core) & (df['section'] == section), 'bot_depth'] += top_depth
        df.loc[(df['core'] == core) & (df['section'] == section), 'ave_depth'] += top_depth
        df.loc[(df['core'] == core) & (df['section'] == section), 'y_m'] += y/1000
        df.loc[(df['core'] == core) & (df['section'] == section), 'x_m'] += x/1000


# %%
# wrap up - reordering columns, saving to csv

# reorder columns
columns_to_move = ['core','section','cut','ave_depth','ave_depth_adj','y_m','x_m']
columns = df.columns.tolist()
for col in columns_to_move:
    columns.remove(col)
columns = columns_to_move + columns
df = df[columns]

# save to csv
df.to_csv(path_to_data + 'master_icpms.csv',index=False)