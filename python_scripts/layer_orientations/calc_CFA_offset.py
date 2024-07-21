#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 08:52:05 2024

Calc offset to CFA sticks

@author: Liam
"""


#%% Import Packages

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
from matplotlib import lines

# weighted percentile statistics
from weighted_percentile import weighted_percentile

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

#%% compute true angles, plot

def truedip(core,side,o_lost,title):

    df = pd.read_pickle('../../data/angles/'+core+'_angles.df')
    for n in ['AC-true-angles','AC-true-scores',
       'DC-true-angles','DC-true-scores',
       'AC-true-orientations','DC-true-orientations']:
        df[n]=None
    
    depth = df['depth'].to_numpy()
    
    
    # loop through all sections and compute drue dip
    for index, row in df.iterrows():
        
        print("Running row "+str(row['section']))
        
        for ACorDC in ['AC']:
            depth = row['depth']
            top_angle = row[ACorDC+'-t-angles']
            top_score = row[ACorDC+'-t-scores']
            side_angle = row[ACorDC+'-'+side+'-angles']
            side_score = row[ACorDC+'-'+side+'-scores']
            top_length = row[ACorDC+'-t-length']
            side_length = row[ACorDC+'-'+side+'-length']
    
    # make plot
    fig,axs = plt.subplots(1,1,figsize=(9,7),dpi=250)
    
    for name,horoff,vertoff in zip(['Dart','OSU1','OSU2'],[1.6,3.2+2,3.2+4+2],[3.2+2]):
        axs.plot(depth)
    
#%% Run

truedip('alhic2302','l',[19.73,20.39,24.28,30.72,39.53,40.47,41.39,42.7],'ALHIC2302')
