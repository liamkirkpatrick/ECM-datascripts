#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:15:03 2024

Angle Calculations on ALHIC2201

@author: Liam
"""

#%% Import packages

# general
import numpy as np
import pandas as pd

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
window = 10

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/first_plot/'
metadata_file = 'metadata.csv'

# ONE BIG PLOT
onebigplot = True
indplots = False

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

sec = set(sections)