#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 15:52:47 2024

Script to:
    1) Calculate true dip and orientation
    2) Make a plot

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

# math
from statsmodels.stats.weightstats import DescrStatsW

#%% Box and Wisker Plot

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


    
#%% compute true angles, plot

def truedip(core,side,o_lost,title):
    # core='alhic2302'
    # side = 'l'
    # o_lost = [19.73,20.39,24.28,30.72,39.53,40.47,41.39,42.7]
    
    df = pd.read_pickle('../../data/angles/'+core+'_angles.df')
    for n in ['AC-true-angles','AC-true-scores',
       'DC-true-angles','DC-true-scores',
       'AC-true-orientations','DC-true-orientations']:
        df[n]=None
    

    
    # loop through all sections and compute drue dip
    for index, row in df.iterrows():
        
        #print("Running row "+str(row['section']))
        
        for ACorDC in ['AC','DC']:
            
            top_angle = row[ACorDC+'-t-angles']
            top_score = row[ACorDC+'-t-scores']
            side_angle = row[ACorDC+'-'+side+'-angles']
            side_score = row[ACorDC+'-'+side+'-scores']
            top_length = row[ACorDC+'-t-length']
            side_length = row[ACorDC+'-'+side+'-length']
            
            dip=[]
            orientation=[]
            score=[]
            for i in range(len(top_angle)):
                a1 = top_angle[i]
                for j in range(len(side_angle)):
                    a2=side_angle[j]
                    score.append((top_score[i]**2)*(side_score[j]**2)*top_length[i]*side_length[j])
                    eps = np.arctan( 1/np.tan(a1*np.pi/180) * np.tan(a2 * np.pi/180)) * 180/np.pi
                    true = np.arctan(np.tan(a1 * np.pi/180) / np.cos(eps*np.pi/180))* 180/np.pi
                    eps = eps+90
                    if true<0:
                        eps = eps+180
                        true = true * -1
                    orientation.append(eps)
                    dip.append(true)
    
            # save to dataframe
            df.at[index,ACorDC+'-true-angles'] = dip
            df.at[index,ACorDC+'-true-scores'] = score
            df.at[index,ACorDC+'-true-orientations'] = orientation
    
    # save dataframe
    df.to_pickle('../../data/angles/'+core+'_angles.df')
    
    # make plot
    fig,axs = plt.subplots(1,2,figsize=(9,7),dpi=250)
    
    # track stats
    central_estimate = []
    
    # loop through and calcualte percentiles
    for index,row in df.iterrows():
        #print("Running row "+str(row['section']))
        for ACorDC,c in zip(['AC'],['r.']):
            dip = np.array(row[ACorDC+'-true-angles'])
            orientation = np.array(row[ACorDC+'-true-orientations'])
            scores = np.array(row[ACorDC+'-true-scores'])
            depth = row['depth']
            
            percentiles = [10,25,50,75,90]
            
            # plot dip
            dip_stats = weighted_percentile(dip, percentiles, weights=scores, interpolation='step')
            wiskerplot(depth,dip_stats,axs[0],ACorDC)
            
            central_estimate.append(dip_stats[2])
            
            # option to plot 
            for j in range(len(scores)):
                #axs[1].plot(orientation[j],depth,c,markersize=(scores[j])*10)
                axs[1].plot(orientation[j],depth,'k.',markersize=2)
                axs[0].plot(dip[j],depth,'k.',markersize=2)

            # alternative approach to orientation stats
            if False:
                orientation_stats = weighted_percentile(orientation, percentiles, weights=scores, interpolation='step')
                wiskerplot(depth,orientation_stats,axs[1],ACorDC)
                
            else:
                sort_index = np.argsort(orientation)
                orientation = orientation[sort_index]
                scores = scores[sort_index]
                spc_threshold = 180
                max_spc = orientation[0] + 360 - orientation[-1]
                bound1 = orientation[-1]
                bound2 = orientation[0]
                for i in range(len(orientation[:-2])):
                    if orientation[i+1]-orientation[i] > max_spc:
                        max_spc = orientation[i+1]-orientation[i]
                        bound1 = orientation[i]
                        bound2 = orientation[i+1]
                if max_spc > spc_threshold:
                    if bound1<bound2:
                        oshift = (orientation + 360-(bound1+bound2)/2)%360
                    else:
                        oshift = orientation
                    orientation_stats = weighted_percentile(oshift, percentiles, weights=scores, interpolation='step')
                    if bound1<bound2:
                        orientation_stats = (orientation_stats+(bound1+bound2)/2)%360
                    
                    wiskerplot(depth,orientation_stats,axs[1],ACorDC)
                #else:
                    
                #    print("    Too small of a gap!")
            

            
    # housekeeping
    axs[0].set_xlabel('Angle (degrees)')
    axs[1].set_xlabel('Orientation (degrees)')
    axs[0].set_ylabel('Depth (m)')
    axs[1].set_ylabel('Depth (m)')
    axs[0].set_title('Dip Angle (Relative to Horizontal)')
    axs[1].set_title('Dip Direction (Relative to Orientation Line)')
    axs[0].grid()
    axs[1].grid()
    axs[0].set_xlim([0,90])
    axs[1].set_xlim([0,360])
    topcushion = (max(df['depth'])-min(df['depth']))/3.5
    axs[0].set_ylim([max(df['depth'])+topcushion/5,min(df['depth'])-topcushion])
    axs[1].set_ylim([max(df['depth'])+topcushion/5,min(df['depth'])-topcushion])
    axs[1].yaxis.tick_right()
    axs[1].yaxis.set_label_position("right")

    vertical_line = lines.Line2D([], [], color='k', marker='|',
                                 linestyle='None', markersize=10,
                                 markeredgewidth=1.5, label='Median')
    # make plot for axis labels
    for a in axs:
        a.plot([-10,-10],[0,1],'k-',label='Median') 
        a.add_patch(Rectangle((-20,0),10,10,
                                facecolor='r',
                                edgecolor='k',
                                linewidth=1.5,
                                label='Weighted Interquartile Range'))
        a.plot([-20,-10],[0,0],'k-',label='Weighted 10%-90% Percentile Range') 
        a.plot([-20],[0],'k.',label='Outlying Estimates') 
    
    # Orientation lost
    if len(o_lost)>0:
        for d in o_lost:
            axs[1].plot([0,360],[d,d],'g-',linewidth=2,label='Orientation Lost')
    

    # make the legend 
    handles, labels = axs[0].get_legend_handles_labels()
    handles.insert(0,vertical_line)
    labels.insert(0,'Weighted Median')
    by_label = dict(zip(labels, handles))
    axs[0].legend(by_label.values(), by_label.keys(),loc='upper right')
    handles, labels = axs[1].get_legend_handles_labels()
    handles.insert(0,vertical_line)
    labels.insert(0,'Weighted Median')
    by_label = dict(zip(labels, handles))
    axs[1].legend(by_label.values(), by_label.keys(),loc='upper right')
    fig.suptitle(title)
    plt.tight_layout()
    fig.savefig('../../../figures/orientations/'+core+'_angleplot.png')
    
    # print stats:
    print(central_estimate)
    est = np.array(central_estimate)
    print("Core: "+core)
    print("    min = "+str(np.min(est)))
    print("    max = "+str(np.max(est)))
    print("    mean = "+str(np.mean(est)))
    print("    std = "+str(np.std(est)))
    
    #%% Make another plot
    
    fig2,axs2 = plt.subplots(1,1,dpi=200)

    weighted_mean=[]
    weighted_std = []
    
    # loop through and calcualte percentiles
    for index,row in df.iterrows():
        #print("Running row "+str(row['section']))
        for ACorDC,c in zip(['AC'],['r.']):
            dip = np.array(row[ACorDC+'-true-angles'])
            orientation = np.array(row[ACorDC+'-true-orientations'])
            scores = np.array(row[ACorDC+'-true-scores'])
            depth = row['depth']
    
                
            weighted_stats = DescrStatsW(dip, weights=scores, ddof=0)
            weighted_mean.append(weighted_stats.mean)
            weighted_std.append(weighted_stats.std)
            
    axs2.plot(abs(weighted_mean-np.mean(weighted_mean)),weighted_std,'k.')
            
    axs2.set_xlabel('Difference from Average Weighted Mean')
    axs2.set_ylabel('Weigthed Standard Deviation')
    axs2.set_title(core)
            
    fig2.savefig('../../../figures/orientations/'+core+'_dipspread.png')
    
    #%% Make another plot
    
    fig3,axs3 = plt.subplots(1,1,dpi=200)

    weighted_mean=[]
    weighted_std = []
    spread = []
    
    # loop through and calcualte percentiles
    for index,row in df.iterrows():
        #print("Running row "+str(row['section']))
        for ACorDC,c in zip(['AC'],['r.']):
            dip = np.array(row[ACorDC+'-true-angles'])
            orientation = np.array(row[ACorDC+'-true-orientations'])
            scores = np.array(row[ACorDC+'-true-scores'])
            depth = row['depth']
    
                
            weighted_stats = DescrStatsW(dip, weights=scores, ddof=0)
            weighted_mean.append(weighted_stats.mean)
            weighted_std.append(weighted_stats.std)
            spread.append(row['10-90 percentile spread'])
            
            
    axs3.plot(spread,weighted_std,'k.')
            
    axs3.set_xlabel('10-90 Percentile Ratio')
    axs3.set_ylabel('Weighted Standard Deviation')
    axs3.set_title(core)
    
            
    fig3.savefig('../../../figures/orientations/'+core+'_dip_percentile.png')
    
    
#%% Run

truedip('alhic2302','l',[10.37,19.73,20.39,24.28,30.72,39.53,40.47,41.39,42.7],'ALHIC2302')
truedip('alhic2201','r',[17.62,17.65,24.08],'ALHIC2201')
