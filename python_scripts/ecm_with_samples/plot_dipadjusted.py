#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:39:20 2024

Plot all data from ALHCI2302 BID full rounds (quarter core measurements)

@author: Liam
"""

#%%
# Import packages 

# general
import numpy as np
import pandas as pd
import os
import glob
import re

# plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# annimation - try 2
#import cv2

# annimation
import moviepy.video.io.ImageSequenceClip
import moviepy.editor as mp
#from PIL import Image
import PIL
PIL.Image.ANTIALIAS = PIL.Image.LANCZOS

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

#%% Function to ensure even dimensions on figure

def ensure_even_dimensions(clip):
    width, height = clip.size
    width = width if width % 2 == 0 else width - 1
    height = height if height % 2 == 0 else height - 1
    return clip.resize(newsize=(width, height))

#%% 
# User Inputs

# set which plots to make
make_singles = False
make_annimations = True

# smoothing window, mm
window = 20

# paths
path_to_data = '../../data/'
path_to_raw = '/Users/Liam/Desktop/UW/ECM/raw_data/'
path_to_figures = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/dip_adjusted_BID_quarter/'
path_to_annimations = '/Users/Liam/Desktop/UW/ECM/2024_structure/figures/deep_annimations/'
metadata_file = 'metadata.csv'
path_to_angle = '../../data/angles/'
path_to_proxies = '../../data/sampling/'


# ONE BIG PLOT
indplots = True

# set sections to run
to_run = ['155_2','156_2','158','159_3']
#to_run = ['159_3']

#%% 
# Read in metadata and import data
meta = pd.read_csv(path_to_data+metadata_file)

#%% 
# import each script as an ECM class item

data = []
cores = []
sections = []
faces = []
ACorDCs = []
for index,row in meta.iterrows():
    
    core = row['core']
    
    if core == 'alhic2302':
        
        section = row['section']
        face = row['face']
        ACorDC = row['ACorDC']
        
        data_item = ECM(core,section,face,ACorDC)
        
        if section in to_run and ACorDC == 'AC' and (face == 'r' or face == 'tr'):
            print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
            
            data_item.rem_ends(10)
            data_item.smooth(window)
            data.append(data_item)
            data_item.norm_all()
            
            cores.append(core)
            sections.append(section)
            faces.append(face)
            ACorDCs.append(ACorDC)

sec = set(sections)

#%% 
# Load angles

core = 'alhic2302'
df = pd.read_csv(path_to_angle+core+'_deepangles_means.csv',index_col=0)

#%% 
# Load Proxies

proxy_list = ['water_iso']
proxy_df = pd.read_csv(path_to_proxies+proxy_list[0]+'/master_'+proxy_list[0]+'.csv',index_col=0)
if len(proxy_list)>0:
    for d in proxy_list[1:]:
        new_df = pd.read_csv(path_to_proxies+d+'/master'+d+'.csv',index_col=0)
        proxy_df = pd.concat([proxy_df, df2], axis=0, join='outer', ignore_index=True)



#%% 
# define plotting function
def plotquarter(y_vec,ycor,d,meas,button,axs,rescale,angle,face,res):
    
    # calculate track width (for plotting)
    width = y_vec[1] - y_vec[0]

    for y in y_vec:
        
        # Pull out data for this track
        idx = ycor==y
        tmeas = meas[idx]
        tbut = button[idx]
        tycor = ycor[idx]
        td = d[idx]

        # downsample ECM to save plotting time (as needed)
        if res != 0:
            int_lo = round(min(td),2)
            int_hi = round(max(td),2)
            depth_interp = np.linspace(int_lo,int_hi,int((int_hi-int_lo)/res)+1)
            meas_interp = np.interp(depth_interp,np.flip(td),np.flip(tmeas))
            but_interp = np.interp(depth_interp,np.flip(td),np.flip(tbut))
            td = depth_interp
            tmeas = meas_interp
            tbut = np.round(but_interp)

        
        if face == 'tr':
            offset = y
        else:
            offset = -y
        cor = offset/1000 * np.tan(angle*np.pi/180)
        td = td + cor
        
        for i in range(len(tmeas)-1):
            
            if tbut[i] == 0:
                axs.add_patch(Rectangle((y-(width-0.2)/2,td[i]),(width-0.2),td[i+1]-td[i],facecolor=my_cmap(rescale(tmeas[i]))))
            # else:
            #     axs.add_patch(Rectangle((y-(width-0.2)/2,td[i]),(width-0.2),td[i+1]-td[i],facecolor='k'))
            
    return()

#%% 
# define function to find unique elements in list
def unique(list1):
 
    # initialize a null list
    unique_list = []
 
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    
    return(unique_list)

#%% 
# Make colormaps

# make colormap - ECM
my_cmap = matplotlib.colormaps['Spectral']

# make colormap
cmap_line = matplotlib.colormaps['tab10']



        



#%%
# Define function

def plot_script(sec,tr_a,r_a,data,proxies,res):
    
    # set data to empty
    AC_t = None
    AC_r = None

    #loop through data 
    for d in data:
        
        # find faces
        if d.core=='alhic2302':
            if d.section==sec:
                if d.face == 'tr':
                    AC_t = d
                if d.face == 'r':
                    AC_r = d                    
    
    # find depth max and minimum
    minvec = []
    maxvec = []
    AC_all = []
    for data_face in [AC_t,AC_r]:
        if data_face != None:
            minvec.append(min(data_face.depth))
            maxvec.append(max(data_face.depth))
            AC_all.extend(data_face.meas)

    # set color scaling
    ACpltmin = np.percentile(AC_all,5)
    ACpltmax = np.percentile(AC_all,95)
    ACrescale = lambda k: (k-ACpltmin) /  (ACpltmax-ACpltmin)
    dmin = min(minvec)-0.1
    dmax = max(maxvec)+0.1

    # Make Figure - ECM

    # Actually make figure (sized for number of proxies to plot)
    ratios = [3,3]
    numplots=2
    if len(proxies)>0:
        ratios.append(1.5)
        numplots=3+len(proxies)
    for i in range(len(proxies)):
        ratios.append(2.5)
    fig, ax = plt.subplots(1, numplots,figsize=(8+len(proxies)*3,8),gridspec_kw={'width_ratios': ratios},dpi=100)
    
    # ECM subplot admin
    ax[0].set_xlim([110, 0])
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_xlim([0,110])
    for a in [ax[0],ax[1]]:
        a.set_ylabel('Depth (m)')
        a.set_xlabel('Distance From Center (mm)',fontsize=6)
        a.set_ylim([dmax, dmin])
    
    # blank subplot admin
    if len(proxies)>0:
        ax[2].axis('off')
        
    # proxy subplot admin
    if len(proxies)>0:
        
        ax[3].set_ylabel('Depth (m)')
        ax[3].yaxis.set_label_position('left')
        ax[3].set_ylim([dmax, dmin])
        
        ax[len(proxies)+2].yaxis.set_label_position('right')
        ax[len(proxies)+2].yaxis.tick_right()
        
        # add axis labels/y-axis limits for all subpots
        for i in range(len(proxies)):
            ax[i+3].set_xlabel(proxies[i])
            ax[i+3].set_title(proxies[i])
            ax[i+3].set_ylim([dmax, dmin])
        
        # remove y-axis labels for middle plots if there are more than 2
        if len(proxies)>2:
            for i in range(len(proxies)-2):
                ax[i+4].yaxis.set_ticks([])
                ax[i+4].set_yticklabels([]) 
    
    # Plot ECM
    for a,data_face,angle in zip([ax[0],ax[1]],[AC_r,AC_t],[r_a,tr_a]):
        if data_face != None:
            if data_face.face == 'r':
                yall = data_face.y_s - data_face.y_left
                yvec = data_face.y_vec - data_face.y_left
            else:
                yall = data_face.y_right - data_face.y_s
                yvec = data_face.y_right - data_face.y_vec
            plotquarter(yvec,
                        yall,
                        data_face.depth_s,
                        data_face.meas_s,
                        data_face.button_s,
                        a,
                        ACrescale,
                        angle,
                        data_face.face,
                        res)
        
    # Plot Samples
    # filter for the correct section
    sec_proxy_df = proxy_df[(proxy_df['section'] == sec)&(proxy_df['core'] == core)]
    # loop through all of the proxy
    prox_cnt = 0
    for proxy in proxies:
        # filter df for this proxy
        df_prox = sec_proxy_df[sec_proxy_df[proxy].notna()]

        # now we have to loop through each stick
        cuts = sec_proxy_df['cut'].unique()
        c_cnt = 0
        for c in cuts:

            # filter for this cut
            df_cut = df_prox[df_prox['cut']==c]

            if len(df_cut)>0:
            
                # pull out the desired info
                top_depth = df_cut['top_depth'].to_numpy()
                bot_depth = df_cut['bottom_depth'].to_numpy()
                ave_depth = (top_depth + bot_depth)/2
                val = df_cut[proxy].to_numpy()
                y = df_cut['y_m'].to_numpy()
                x = df_cut['x_m'].to_numpy()
                # need to impliment something to check the cut is consistent through the stick 
                # - here I assume the first sample is representative. Should be true regardless
                # but is still worth checking
                y_shift = -1 * y[0] * np.tan(r_a*np.pi/180)
                x_shift = 1 * x[0] * np.tan(tr_a*np.pi/180)
                shift = y_shift + x_shift

                # plot on dedicated subplot
                ax[3+prox_cnt].plot(val,ave_depth+shift,color=cmap_line(c_cnt),label=c)
                for i in range(len(top_depth)):

                    ax[3+prox_cnt].plot([val[i],val[i]],[top_depth[i]+shift,bot_depth[i]+shift],color=cmap_line(c_cnt))

                ax[3+prox_cnt].legend()


            # increment Counter
            c_cnt+=1


        # increment counter
        prox_cnt+=1


    # Plot housekeeping
    fig.suptitle('ALHIC2302 - '+sec+' - tr='+str(round(tr_a))+'- r='+str(round(r_a)))
    ax[0].set_title('AC - Right')
    ax[1].set_title('AC - Top')
    for a in ax:
        a.grid(True)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0)

    # add colorbar (width fixed to first two subplots)
    ACcbar_ax = fig.add_axes([ax[0].get_position().x0,-0.05,ax[1].get_position().x1-ax[0].get_position().x0,0.05])
    ACnorm = matplotlib.colors.Normalize(vmin=ACpltmin,vmax=ACpltmax)
    ACcbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=ACnorm, cmap=my_cmap),cax=ACcbar_ax,
                    orientation='horizontal',label='Current (normalized)')

    return fig

#%% 
# Loop through and plot (just endpoints)

if make_singles:

    # print update
    print("Running Single Plots")

    # Specify proxies to plot
    proxies=['d18O','dxs']

    # set resolution (0 is default of measurment resolution)
    res = 0

    # loop through sections
    for sec in unique(sections):

        # set angles to run
        tr_angle = [0,0,float(df.at[sec,'AC-tr-mean-angle']),float(df.at[sec,'AC-tr-mean-angle'])]
        r_angle = [0,float(df.at[sec,'AC-r-mean-angle']),0,float(df.at[sec,'AC-r-mean-angle'])]
        
        for tr_a, r_a in zip(tr_angle,r_angle):

            # print update
            print("    Running Section "+sec+' - tr='+str(round(tr_a))+' - r='+str(round(r_a)))
            
            # run plot script
            fig = plot_script(sec,tr_a,r_a,data,proxies,res)
            
            #save figure
            fname = path_to_figures +'alhic2302-'+sec+'_tr='+str(round(tr_a))+'_r='+str(round(r_a))+'.png'
            fig.savefig(fname,bbox_inches='tight')
            plt.close(fig)

#%% 
# Make some annimations

if make_annimations:

    # print update
    print("Plotting Annimations")

    # Specify proxies to plot
    proxies=['d18O','dxs']

    # set resolution
    res = 0.002

        # loop through sections
    for sec in unique(sections):

        # set folder path
        folder_path = path_to_annimations+'/'+sec+'/'

        # print update
        print("    Running Section "+sec)

        # delete files in folder
        png_files = glob.glob(os.path.join(folder_path, "*.png"))
        # Loop through the list of .png files and delete each one
        for file_path in png_files:
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")

        # set counter
        fig_cnt = 1
        
        # make vectors of tr angles and r angles
        angle_res = 0.2
        tr_final = float(df.at[sec,'AC-tr-mean-angle'])
        r_final = float(df.at[sec,'AC-r-mean-angle'])
        phase1_cnt = round(abs(tr_final/angle_res))+1
        phase2_cnt = round(abs(r_final/angle_res))+1
        tr_angle = np.concatenate((np.zeros(phase1_cnt),np.linspace(0,tr_final,phase2_cnt)))
        r_angle = np.concatenate((np.linspace(0,r_final,phase1_cnt),np.ones(phase2_cnt)+r_final))
        
        #tr_angle = [0,0,float(df.at[sec,'AC-tr-mean-angle']),float(df.at[sec,'AC-tr-mean-angle'])]
        #r_angle = [0,float(df.at[sec,'AC-r-mean-angle']),0,float(df.at[sec,'AC-r-mean-angle'])]
        
        # loop through angles
        for tr_a, r_a in zip(tr_angle,r_angle):

            # print update
            print("        Plotting subplot "+str(fig_cnt)+" of "+str(len(tr_angle)))
            
            # run plot script
            fig = plot_script(sec,tr_a,r_a,data,proxies,res)
            
            #save figure
            fname = path_to_annimations+'/'+sec+'/alhic2302-'+sec+'-'+str(fig_cnt)+'.png'
            fig.savefig(fname,bbox_inches='tight')
            plt.close(fig)

            # increment counter
            fig_cnt += 1

        # make movie
 
        def get_png_files(folder_path):
            return [file for file in os.listdir(folder_path) if file.endswith('.png')]

        image_names = get_png_files(folder_path)
        image_names = sorted(image_names, key=lambda x: int(re.findall(r'\d+', x)[-1]))
        image_files = [os.path.join(folder_path,img)
                    for img in image_names]
    

        clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files,fps=30)



        clip = ensure_even_dimensions(clip)
        clip.write_videofile(folder_path+sec+'_movie.mp4',fps=24)

# %%
