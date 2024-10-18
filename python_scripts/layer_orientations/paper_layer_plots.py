#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 8, 2024

Script to:
    1. Compute the true dip and orientation of the core
    2. Plot the results (for ECM paper)

@author: Liam
"""

#%% 
# Import Packages

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

# math
from statsmodels.stats.weightstats import DescrStatsW

#%% Load Data
# Note - must run get_truedip first to actually calculate the true dip

alhic2302 = pd.read_pickle('../../data/angles/alhic2302_angles.df')
alhic2201 = pd.read_pickle('../../data/angles/alhic2201_angles.df')

#%% 
# calculate true dip percentiles for alhic2302 and alhic2201

# define function for weighted_percentile
def weighted_percentile(values, weights, percentiles):
    
    # Convert percentiles to fractions
    percentiles = [p / 100 for p in percentiles]

    # Convert values and weights to numpy arrays
    values = np.array(values)
    weights = np.array(weights)
    percentiles = np.array(percentiles)

    # Check for NaN values and remove them
    nan_mask = np.isnan(values) | np.isnan(weights)
    if np.any(nan_mask):
        print("NaN values found and will be removed.")
    values = values[~nan_mask]
    weights = weights[~nan_mask]
    
    # Sort values and weights by values
    sorted_indices = np.argsort(values)
    sorted_values = values[sorted_indices]
    sorted_weights = weights[sorted_indices]
    
    # Compute the cumulative sum of weights
    cumulative_weights = np.cumsum(sorted_weights)
    
    # Normalize the cumulative weights to get the cumulative distribution
    cumulative_distribution = cumulative_weights / cumulative_weights[-1]
    
    # Interpolate to find the weighted percentiles
    weighted_percentiles = np.interp(percentiles, cumulative_distribution, sorted_values)
    
    return weighted_percentiles

# define function to calculate percentiles from dataframe with true dip angles
def calc_percentiles(df,percentiles):

    # initialize list to store dip statistics
    dip_stats = []
    # define if we are looking at AC or DC (AC is default)
    ACorDC = 'AC'

    # loop through each row in the dataframe
    for index,row in df.iterrows():

        # pull out approriate data
        dip = np.array(row[ACorDC+'-true-angles'])
        scores = np.array(row[ACorDC+'-true-scores'])
        
        # plot dip
        a = weighted_percentile(dip, scores, percentiles)
       # weighted_percentile(dip, percentiles, weights=scores, interpolation='step')
        dip_stats.append(a)

        # check
        for i in range(len(a)-1):
            if a[i]>a[i+1]:
                print('Error: Percentiles are not in order in section '+row['section'])

    df['dip_percentiles']  = dip_stats
    return df

#deine percentiles
percentiles = [10,25,50,75,90]

# calculate true dip percentiles for alhic2201 and alhic2302
alhic2201 = calc_percentiles(alhic2201,percentiles)
alhic2302 = calc_percentiles(alhic2302,percentiles)

#%% 
# Make Plot of depth vs dip

def wiskerplot(d,q,axs,ACorDC):
    
    h = 0.8
    linewidth=3
    
    if ACorDC == 'AC':
        c = 'k'
        o = 0
    else:
        c = 'b'
        o = 0
        
    #box
    if q[1]<q[3]:
        axs.add_patch(Rectangle((q[1],d-h/2+o),(q[3]-q[1]),
                                h,
                                facecolor='r',
                                edgecolor='k',
                                linewidth=1))
    else:
        axs.add_patch(Rectangle((q[3],d-h/2+o),(360-q[3]),
                                h,
                                facecolor='r',
                                edgecolor='k',
                                linewidth=1))
        axs.add_patch(Rectangle((0,d-h/2+o),(q[1]),
                                h,
                                facecolor='r',
                                edgecolor='k',
                                linewidth=1))
    
    # centerline
    if q[0]<q[4]:
        axs.plot([q[0],q[-1]],[d+o,d+o],c,linewidth=linewidth)
    else:
        axs.plot([q[0],360],[d+o,d+o],c,linewidth=linewidth)
        axs.plot([0,q[4]],[d+o,d+o],c,linewidth=linewidth)

    # plot vertical lines at each quadrant
    axs.plot([q[0],q[0]],[d-h/4+o,d+h/4+o],c,linewidth=linewidth)
    #axs.plot([q[1],q[1]],[d-h/2+o,d+h/2+o],c,linewidth=linewidth)
    axs.plot([q[2],q[2]],[d-h/2+o,d+h/2+o],c,linewidth=linewidth)
    #axs.plot([q[3],q[3]],[d-h/2+o,d+h/2+o],c,linewidth=linewidth)
    axs.plot([q[4],q[4]],[d-h/4+o,d+h/4+o],c,linewidth=linewidth)


# make figure
fig, ax = plt.subplots(1,2,figsize=(8,10),dpi=300)

# axis labels
for a in ax:
    a.set_ylabel('Depth (m)')
    a.set_xlabel('Dip Angle (degrees)')
    a.set_xlim(0,90)
    a.grid(True)
ax[0].set_title('ALHIC2201')
ax[1].set_title('ALHIC2302')
ax[0].invert_yaxis()
ax[1].invert_yaxis()
ax[0].set_ylim(23,8)
ax[1].set_ylim(47,8)

# plot data
for index,row in alhic2201.iterrows():

    d = row['depth']
    dip = np.array(row['AC-true-angles'])
    orientation = np.array(row['AC-true-orientations'])
    scores = np.array(row['AC-true-scores'])
    wiskerplot(d,row['dip_percentiles'],ax[0],'AC')
for index,row in alhic2302.iterrows():
    d = row['depth']
    dip = np.array(row['AC-true-angles'])
    orientation = np.array(row['AC-true-orientations'])
    scores = np.array(row['AC-true-scores'])
    wiskerplot(d,row['dip_percentiles'],ax[1],'AC')


#%%
# Make a plot showing dip direction

# function to calculate circular mean
def calc_circular_mean(df):

    for index,row in df.iterrows():
        
        # pull out approriate data
        angles = np.array(row['AC-true-orientations'])
        weights = np.array(row['AC-true-scores'])
        
        # convert to radians
        angles = np.radians(angles)
        
        # compute the mean
        x = np.sum(np.cos(angles)*weights)
        y = np.sum(np.sin(angles)*weights)
        
        # convert back to degrees
        mean = np.degrees(np.arctan2(y,x))

        # convert to 0-360
        if mean < 0:
            mean = mean + 360
        
        # store the circular mean
        df.loc[index,'orientation'] = mean
    
    return df

alhic2201 = calc_circular_mean(alhic2201)
alhic2302 = calc_circular_mean(alhic2302)

# define core breaks
alhic2201_breaks = [17.62,17.65,24.08]
alhic2302_breaks = [10.37,19.73,20.39,24.28,30.72,39.53,40.47,41.39,42.7]

# make figure
fig, ax = plt.subplots(1,1,figsize=(4,6),dpi=300)

ax.set_ylabel('Depth (m)')
ax.set_xlabel('Dip Direction (degrees)')
ax.grid(True)
ax.set_title('ALHIC2302 Dip Direction')
ax.set_xlim(0,360)
ax.set_ylim(47,10)

# plot dip directions
cnt = 0
for index,row in alhic2302.iterrows():
    d = row['depth']
    if cnt==0:
        ax.plot(row['orientation'],d,'ro',markersize=3,label='Dip Direction')
    else:
        ax.plot(row['orientation'],d,'ro',markersize=3)
    cnt+=1

# plot core breaks
cnt = 0
for b in alhic2302_breaks:
    if cnt ==0:
        ax.plot([0,360],[b,b],'k--',label='Orientation Lost')
    else:
        ax.plot([0,360],[b,b],'k--')
    cnt+=1

# add legend
ax.legend()


#%%
# Make a plot showing t-test of trend

def t_test_trend(df):

    # pull out median dip
    dip = np.array(df['dip_percentiles'].apply(lambda x: x[2]))

    # pull out depth
    depth = np.array(df['depth'])

    # perform t-test
    slope, intercept, r_value, p_value, std_err = stats.linregress(depth,dip)

    return slope,intercept,p_value,depth,dip

# calculate t-test for alhic2201 and alhic2302
slope2201,intercept2201,p_value2201,depth2201,dip2201 = t_test_trend(alhic2201)
slope2302,intercept2302,p_value2302,depth2302,dip2302 = t_test_trend(alhic2302)

# make figure
fig, ax = plt.subplots(2,1,figsize=(7,7),dpi=300)
for a in ax:
    a.set_xlabel('Depth (m)')
    a.set_ylabel('Dip Angle (degrees)')
    a.set_xlim(0,90)
    a.grid(True)
ax[0].set_title('ALHIC2201')
ax[1].set_title('ALHIC2302')
ax[0].set_xlim(8,23)
ax[1].set_xlim(8,47)

# plot data
for ax,slope,intercept,p_value,depth,dip in zip(ax,[slope2201,slope2302],[intercept2201,intercept2302],[p_value2201,p_value2302],[depth2201,depth2302],[dip2201,dip2302]):

    # plot data
    ax.plot(depth,dip,'ro',markersize=3,label='Median Dip Angle')

    # plot trendline
    label = f'Trendline (Slope: {slope:.2f}, \n Pval = {p_value:.3f})'
    ax.plot([0,max(depth)+3],[intercept,intercept+slope*(max(depth)+3)],'k--', label=label)

    # add p-value
    #ax.text(0.1,0.9,'p = '+str(round(p_value,3)),transform=ax.transAxes)

    ax.legend()
plt.tight_layout()

# %%
# save again