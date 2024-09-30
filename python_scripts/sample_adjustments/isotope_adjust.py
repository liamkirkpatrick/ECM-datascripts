#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 17:37:45 2024

@author: Liam
"""

#%% load packages

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%% Set filepaths

path_to_data = '../../data/sampling/coordinates/'
warning_threshold = [0.5,1.5]

#%% Load data

df = pd.read_excel(path_to_data+'isotope_depths.xlsx', sheet_name='sample_thicknesses',index_col=0)
meta = pd.read_excel(path_to_data+'isotope_depths.xlsx', sheet_name='metadata',index_col=0)

#%% Function to add column to DF (thanks chatGPT)

def add_column_to_df(df, new_list, column_name):
    # Find the current number of rows in the DataFrame
    current_rows = df.shape[0]
    
    # Find the length of the new list
    new_list_length = len(new_list)
    
    # If the new list is longer, extend the DataFrame with NaN values
    if new_list_length > current_rows:
        df = df.reindex(range(new_list_length))
    
    # If the new list is shorter, pad it with NaN to match the DataFrame size
    elif new_list_length < current_rows:
        new_list = new_list + [np.nan] * (current_rows - new_list_length)
    
    # Add the new list as a new column
    df[column_name] = new_list
    
    return df

#%% Apply corrections

sticks = df.columns
output = pd.DataFrame()

for s in df.columns[0:]:
    
    print("Running stick "+s)
    
    
    # get stick metadata
    top_depth = meta['section_top_depth_m'][s]
    offset = meta['offset_mm'][s]
    length = meta['length'][s]
    
    
    # extract sample lengths and calculate total length
    s_len = np.array(df[df[s].notna()][s])
    total_s_len = sum(s_len)
    print("    Length = "+str(len(s_len)))
    
    # calculate difference between sample length sum and total length
    diff = length - total_s_len
    
    # calculate width of average cut, update reader
    cut_width = diff / (len(s_len)-1)
    if cut_width > max(warning_threshold) or cut_width < min(warning_threshold):
        print("    **** Warning: Cut Exceeds Expected Threshold ****")
    print("    Stick Length   = "+str(length))
    print("    Sample Lengths = "+str(total_s_len))
    print("    Ave Cut Width  = "+str(cut_width))
    
    # calculate top/bottom/mid depths,
    s_td = []
    s_bd = []
    for i in range(len(s_len)):
        if i==0:
            s_td.append(np.round(top_depth + offset/1000,3))
        # elif i == 1:
        #     s_td.append(np.round(top_depth + (offset+s_len[0]+cut_width)/1000,3))
        else:
            s_td.append(np.round(top_depth + (offset+sum(s_len[0:i])+(i-1)*cut_width)/1000,3))
        s_bd.append(np.round(top_depth + (offset+sum(s_len[0:i+1])+i*cut_width)/1000,3))
    
    # add results to output    
    output = add_column_to_df(output, s_td, s+'_td')
    output = add_column_to_df(output, s_bd, s+'_bd')

# write to excel sheet
with pd.ExcelWriter(path_to_data+'isotope_depths.xlsx',
                    mode='a', engine='openpyxl',
                    if_sheet_exists='replace') as writer:
    output.to_excel(writer, sheet_name='Depths')
        
