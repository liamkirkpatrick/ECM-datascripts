#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 22:31:42 2024

@author: Liam
"""

#%% Import packages

import pandas as pd
import numpy as np
import math

#%% User inputs

path_to_data = '../../data/'
metadata = 'metadata.csv'

#%% Define ECM Class

class ECM:
    
    def __init__(self,core,section,face,ACorDC):
        
        # open metadata csv
        meta = pd.read_csv(path_to_data+metadata)
        row = meta.loc[meta['core']==core]
        row = row.loc[row['section']==section]
        row = row.loc[row['face']==face]
        row = row.loc[row['ACorDC']==ACorDC]
        
        
        # assign core components
        self.time = row['time'].values[0]
        self.y_left = row['Y_left'].values[0]
        self.y_right = row['Y_right'].values[0]
        self.core = row['core'].values[0]
        self.section = row['section'].values[0]
        self.face = row['face'].values[0]
        self.ACorDC = row['ACorDC'].values[0]
        
        # open core csv
        fname = self.core+'-'+self.section+'-'+self.face+'-'+self.ACorDC+'.csv'
        raw = pd.read_csv(path_to_data+self.core+'/'+fname)
        print(raw)
        
        # assign vectors
        self.meas = raw['meas'].to_numpy()
        self.y = raw['Y_dimension(mm)'].to_numpy()
        self.button = raw['Button'].to_numpy()
        self.depth = raw['True_depth(m)'].to_numpy()
        self.y_vec = np.unique(self.y)
        
    def smooth(self,window):
        # takes as input smoothing window (in mm)
        
        # convert window to m
        window=window/1000
        
        # get spacing between points
        vec = self.depth[self.y==self.y_vec[0]]
        dist = []
        for i in range(len(vec)-1):
            dist.append(vec[i+1]-vec[i])
        
        # make vector of depths to interpolate onto
        dvecmin = min(self.depth) + 2*window/3
        dvecmax = max(self.depth) - 2*window/3
        dvec_num = math.floor((dvecmax-dvecmin) / abs(np.mean(dist)))
        depth_vec = np.linspace(dvecmin,dvecmax,dvec_num)
        #self.dvec = np.flip(depth_vec)
        
        # make empty smooth vectors
        depth_smooth = []
        meas_smooth = []
        button_smooth = []
        y_smooth = []
        
        # loop through all tracks
        for y in self.y_vec:
            
            # index within track
            idx = self.y == y
            
            dtrack = self.depth[idx]
            mtrack = self.meas[idx]
            btrack = self.button[idx]
            
            # loop through all depths
            for d in depth_vec:
                
                # find index of all points within window of this depth
                didx = (dtrack >= d-window/2) * (dtrack <= d+window/2)
                
                # save values
                depth_smooth.append(d)
                meas_smooth.append(np.mean(mtrack[didx]))
                y_smooth.append(y)
                if sum(btrack[didx])>0:
                    button_smooth.append(1)
                else:
                    button_smooth.append(0)
                
        # save smooth values
        self.depth_s = np.flip(np.array(depth_smooth))
        self.meas_s = np.flip(np.array(meas_smooth))
        self.button_s = np.flip(np.array(button_smooth))
        self.y_s = np.flip(np.array(y_smooth))
        
    def rem_ends(self,clip):
        
        # convert clip to m
        clip = clip/1000
        
        # find index within clip
        dmin = min(self.depth)
        dmax = max(self.depth)
        idx = (self.depth>=dmin+clip) * (self.depth<=dmax-clip)
        
        self.meas = self.meas[idx]
        self.y = self.y[idx]
        self.button = self.button[idx]
        self.depth = self.depth[idx]
        
        
        # check if smooth exists
        if hasattr(self,'y_s'):
            
            # find index within clip
            dmin = min(self.depth_s)
            dmax = max(self.depth_s)
            idx = (self.depth_s>=dmin+clip) * (self.depth_s<=dmax-clip)
            
            self.meas_s = self.meas_s[idx]
            self.y_s = self.y_s[idx]
            self.button_s = self.button_s[idx]
            self.depth_s = self.depth_s[idx]
        
#%% Test

if __name__ == "__main__":
    
    test = ECM('alhic2201','10_1','t','AC')
    
    test.smooth(1)
    
    test.rem_ends(1)
    
    test.smooth(1)
    
    # print(ACtest.core)
    
    # # test spacing between ac points
    # vec = ACtest.depth[ACtest.y==ACtest.y_vec[0]]
    # ACdist = []
    # for i in range(len(vec)-1):
    #     ACdist.append(vec[i+1]-vec[i])
        
    # DCtest = ECM('alhic2201','10_1','t','DC')
    
    # print(ACtest.core)
    
    # # test spacing between ac points
    # vec = DCtest.depth[DCtest.y==DCtest.y_vec[0]]
    # DCdist = []
    # for i in range(len(vec)-1):
    #     DCdist.append(vec[i+1]-vec[i])