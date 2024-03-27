#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:23:35 2024

This script reads new data from march '24 trip to ICF, and saves to master 
file structure as .npy files

@author: Liam
"""

#%% Import packages

import numpy as np
import pandas as pd
import os


#%% User Inputs

path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/first_plot/'
metadata_file = 'metadata.csv'

dates = ['2024-03-19','2024-03-20','2024-03-21','2024-03-22']

# set flags in file and corresponding header in master csv
flag_dict = {'AC Collect Speed: ':'AC_col_sp',
'DC Collect Speed: ':'DC_col_sp',
'DC Voltage: ':'DC_volt',
'Note: ':'note',
'mm per encode step: ':'mm_per_encode_step',
'Number of Expected tracks: ':'num_tracks',
'ACDC offset: ':'ACDC_offset',
'Laser offset: ':'laser_offset',
'Y Left: ':'Y_left',
'Y Right: ':'Y_right',
'AC edgespace ':'AC_edgespace',
'DC edgespace ':'DC_edgespace',
'Index Mark (raw - not laser corrected): ':'idx1_raw',
'Index Mark Relative Depth: ':'idx1_rel',
'Index Mark 2 Relative Depth: ':'idx2_rel',
'Index Mark 3 Relative Depth: ':'idx3_rel',
'(first) Index Mark Absolute Depth: ' : 'idx_abs',
'X min Position (raw - not laser corrected): ':'xmin',
'X max Position (raw - not laser corrected): ': 'xmax'
    }

#
extra_headers = ['time','core','section','face','ACorDC']

#%% Create / read in csv with info

# sort headers / keys
headers = list(flag_dict.values())
flags = flag_dict.keys()
for s in extra_headers:
    headers.append(s)

# check if metadata file already exists. If not, create dataframe
if os.path.exists(path_to_data+metadata_file):
    
    # Read the CSV file into a pandas dataframe
    df = pd.read_csv(path_to_data+metadata_file)
    
else:
    df = pd.DataFrame(columns=headers)

#%% Get list of all file names

txt_files = []
for date in dates:
    folder_path = path_to_raw + date
    # Check if the folder exists
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        # Get the list of all .txt files in the folder and add their paths to the txt_files list
        txt_files.extend([os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.txt')])



#%% Populate dataframe

for f in txt_files:
    
    # open file
    
    
    vals=[]
    # extract values from file
    vals = []
    with open(f, 'r') as file:
        
        cnt = 0
        flags = list(flag_dict.keys())
        while cnt<len(flag_dict.keys()):
            for line in file:
                for flag in flags:
                    if flag in line:
                        flags.remove(flag)
                        vals.append(line[len(flag):-6])
                        cnt+=1
                        
            # on last line in file, check for AC or dc
            lp = line.split(',')
            if lp[3]=='--':
                ACorDC = 'AC'
                print('AC')
            elif lp[4]=='--':
                ACorDC = 'DC'
                print('DC')
            else:
                print("ERROR - not AC or DC")
                ACorDC = 'ERROR'
        
                    
        
    # now add on extra headers not in the flags dict
    path = f.split('/')
    parts = path[-1].split('-')
    vals.append(parts[0]+'-'+parts[1]+'-'+parts[2]+'-'+parts[3]+'-'+parts[4])# time
    vals.append(parts[5]) # core
    vals.append(parts[6]) # section
    vals.append(parts[7][:-4]) # face
    vals.append(ACorDC)
    
    # add to df
    data_dict = dict(zip(headers,vals))
    df = pd.concat([df,pd.DataFrame([data_dict])], ignore_index=True)

#%% Save metadata CSV

# rearange so the columns I want are at the front
h = list(df.columns)
front = ['core','time','section','face','ACorDC']
front.reverse()
for f in front:
    h.remove(f)
    h.insert(0,f)
df = df[h]

# drop duplicates (includes time, so only drops duplicates of the SAME RUN)
df.drop_duplicates(subset=front,keep='last')

# save
df.to_csv(path_to_data+metadata_file)

        