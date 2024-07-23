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
metadata_file = 'metadata.csv'

dates = ['2024-03-19','2024-03-20','2024-03-21','2024-03-22','2024-04-08',
          '2024-04-09','2024-04-10','2024-04-11','2024-04-12','2024-04-15',
          '2024-04-16','2024-04-17','2024-04-18','2024-04-19',"2024-07-09",
          "2024-07-10",
          "2024-07-11",
          "2024-07-12",
          "2024-07-15","2024-07-16","2024-07-17","2024-07-18","2024-07-19",'2024-07-22']

#dates = [']

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
'X max Position (raw - not laser corrected): ': 'xmax'}

#header line
hline = 'Y_dimension(mm),X_dimension(mm),Button,AC,DC,True_depth(m)'

# extra headers for dataframe
extra_headers = ['time','core','section','face','ACorDC','header','filename']

#%% Create / read in csv with info

# sort headers / keys
headers = list(flag_dict.values())
flags = flag_dict.keys()
for s in extra_headers:
    headers.append(s)

#check if metadata file already exists. If not, create dataframe
if os.path.exists(path_to_data+metadata_file):
    
    # Read the CSV file into a pandas dataframe
    df = pd.read_csv(path_to_data+metadata_file)
    
    # rearange so the columns I want are at the front
    df = df[headers]
    
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
    
    
    # extract values from file
    vals = []
    with open(f, 'r') as file:
        
        cnt = 0
        flags = list(flag_dict.keys())
        lcnt = 0
        for index,line in enumerate(file,1):
            for flag in flags:
                if flag in line:
                    flags.remove(flag)
                    vals.append(line[len(flag):-6])
            if hline in line:
                header = index
            
                    
        # on last line in file, check for AC or dc
        lp = line.split(',')
        if lp[4]=='--':
            ACorDC = 'AC'
            print('AC')
        elif lp[3]=='--':
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
    section = parts[6]
    
    if len(parts)==8:
        vals.append(parts[6]) # section
        vals.append(parts[7][:-4]) # face
        
    else:
        vals.append(parts[6][:-4]) # section
        vals.append('half')
    vals.append(ACorDC)
    vals.append(header)
    vals.append(f)
    
    # Now, for alhic2302, look at the 
    
    
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
df = df.drop_duplicates(subset='time')

# sort values
df = df.sort_values(['core','section','face','ACorDC'], ascending=[True,True,True,True])

# save
df.to_csv(path_to_data+metadata_file)

#%% Save CSV with edited depth


# loop through all rows in dataframe 
for index,row in df.iterrows():
    
    
    # read in raw data
    raw = pd.read_csv(row['filename'],header = row['header']-1)
    
    # rename AC/DC to meas
    raw['meas'] = raw[row['ACorDC']]
    raw = raw.drop(['AC', 'DC'], axis=1)

    # convert x to depth
    raw['True_depth(m)'] = float(row['idx_abs']) + (float(row['xmax']) - float(row['idx1_raw']) - raw['X_dimension(mm)']) / 1000
    
    # drop X dimension
    raw = raw.drop(['X_dimension(mm)'],axis=1)
    
    fname = row['core']+'-'+str(row['section'])+'-'+row['face']+'-'+row['ACorDC']+'.csv'
    
    raw.to_csv(path_to_data+row['core']+'/'+fname,index=False)
    
        