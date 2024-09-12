#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 10:15:54 2024

This script will allow me to alter the "button" field in proccessed allan
Hills ECM data files.

It includes a GUI for viewing the data and making adjustments

@author: Liam
"""

#%% Import packages

# basic packages
import pandas as pd
import numpy as np

# plotting
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
from matplotlib.patches import Rectangle
matplotlib.use('TkAgg')
global my_cmap
my_cmap = matplotlib.colormaps['Spectral']

# gui
import PySimpleGUI as sg

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

global fig

#%% User Inputs
path_to_data = '../../data/'

update_file = True

#%% return coorinates when I click on the graph
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata

    global coords
    coords = iy
    
    fig.canvas.mpl_disconnect(cid)

    return iy

#%% Handle Plots

# draw figure within GUI
def draw_figure(canvas, figure):
   tkcanvas = FigureCanvasTkAgg(figure, canvas)
   tkcanvas.draw()
   tkcanvas.get_tk_widget().pack(side='top', fill='both', expand=1)
   return tkcanvas

# delete figure
def delete_figure(figure_agg):
    figure_agg.get_tk_widget().forget()
    plt.close('all')

# make plot
def makeplot(ymin,lmin,lmax,xmin,xmax,d_window,data,currtrack,df):
    
    # set up colorbar
    pltmin = np.percentile(d.meas,5)
    pltmax = np.percentile(d.meas,95)
    rescale = lambda k: (k-pltmin) /  (pltmax-pltmin)


    yspc = data.y_vec[1]-data.y_vec[0]
    cmap = matplotlib.colormaps.get_cmap('coolwarm')
    
    fig, ax = plt.subplots(1, 2,figsize=(10, 7), dpi=100)
    
    # plot axis 1 pcolormesh
    if True:
        # Example vectors
        x = data.y_s
        y = data.depth_s
        z = data.meas_s
        x_unique = np.unique(x)
        y_unique = np.unique(y)
    
        # Create a meshgrid
        X, Y = np.meshgrid(x_unique, y_unique)
    
        # Map z to the grid
        Z = np.full(X.shape, np.nan)  # Initialize with NaNs
        for i in range(len(z)):
            ix = np.where(x_unique == x[i])[0]
            iy = np.where(y_unique == y[i])[0]
            Z[iy, ix] = z[i]
    
        # Plotting using pcolormesh
        ax[0].pcolormesh(X, Y, Z, shading='auto',cmap=my_cmap)
        #ax[0].colorbar()  # Show color scale
    
    t = np.arange(0, 3, .01)
    # plot one line for each track (so each independent y-vector
    for i in range(len(data.y_vec)):
        ind = data.y_s==data.y_vec[i]
        
        if data.y_vec[i]==currtrack:
            linewidth = 7
        else:
            linewidth = 2
        
        ax[1].plot(data.meas_s[ind],data.depth_s[ind],linewidth=linewidth,color=cmap(i/len(data.y_vec)))
        
        # overlay button
        if hasattr(data, 'button_raw'):
            ind_ns = data.y==data.y_vec[i]
            idx_button_raw = ind_ns * (data.button_raw==1)
            ax[1].plot(data.meas[idx_button_raw],data.depth[idx_button_raw],'g.',markersize=linewidth*1.1)
        idx_button_s = ind * data.button_s==1
        ax[1].plot(data.meas_s[idx_button_s],data.depth_s[idx_button_s],'k.',markersize=linewidth)
        
        
    
    # plot lmin lines
    cnt = 0
    for x in lmin:
        ax[1].plot([min(data.meas_s), max(data.meas_s)],[x,x],'k')
        ax[0].plot([currtrack-yspc/2, currtrack+yspc/2],[x,x],'k')
        cnt +=1
        
    # plot lmax line and patch
    cnt = 0
    for x in lmax:
        
        ax[1].plot([min(data.meas_s), max(data.meas_s)],[x,x],'g')
        ax[1].add_patch(Rectangle((min(data.meas_s),lmin[cnt]),max(data.meas_s)-min(data.meas_s),x-lmin[cnt],facecolor=(1, 0, 0, 0.2)))
        ax[0].plot([currtrack-yspc/2, currtrack+yspc/2],[x,x],'g')
        ax[0].add_patch(Rectangle((currtrack-yspc/2,lmin[cnt]),yspc,x-lmin[cnt],facecolor=(1, 0, 0, 0.2)))
        cnt +=1
        
    # plot square around current track
    ax[0].add_patch(Rectangle((currtrack-yspc/2,min(data.depth_s)),
                              yspc,max(data.depth_s)-min(data.depth_s),
                              edgecolor='black',facecolor='none'))

    ax[1].set_xlim(xmin,xmax)
    ax[1].set_ylim(ymin,ymin+d_window)
    ax[1].set_ylabel('Depth (m)')
    ax[1].set_xlabel('Conductivity (amps)')
    ax[1].set_title('Individual Curves')
    
    ax[0].set_ylim(ymin,ymin+d_window)
    ax[0].set_xlim(data.y_left,data.y_right)
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Conductivity (amps)')
    ax[0].set_title('Top View')
    
    fig.tight_layout()
    
    return(fig)

#%% make GUI

def make_gui():
    
    # Create GUI
    first_col = [
       [sg.Canvas(key='-CANVAS-')],
       [sg.RealtimeButton(sg.SYMBOL_DOWN, key='-UP-'),
            sg.Text(' '),sg.RealtimeButton(sg.SYMBOL_UP, key='-DOWN-'),
            sg.Button(button_text='Resize x-axis',key='-RESCALE-'),
            sg.Button(button_text='Zoom in',key = '-+Z-'),
            sg.Button(button_text='Zoom out',key = '--Z-'),
            sg.Button('Show Button',key='-BUTTON-')]
       ]

    h = ['Min Location','Max Location']
    tbl1 = sg.Table(values = [['-','-']], headings = h,key='-TBL-')

    seccond_col = [
        [sg.Text('Current File: '),sg.Text(size=(25, 1), key='-FILE-', pad=(1, 1))],
        [sg.Text('Current Track: '),sg.Text(size=(25, 1), key='-TRACK-', pad=(1, 1))],
        [sg.Button(button_text='Go To Next Track',key = '-NEXTTRACK-')],
        [sg.Button(button_text='Go To Next File (save current work)',key = '-NEXTFILE-')],
        [sg.Button(button_text='Go To Next File (do not save this run)',key = '-NEXTFILE_NOSAVE-')],
        [sg.Quit(focus=True)],
        [sg.Text('_'*30)],
        [tbl1],
        [sg.Button('Delete Last',key = '-DEL-')],
        [sg.Text('_'*30)],
        [sg.Text('Status:'),sg.Text(size=(35,1),key='-STATUS-',justification = 'c')],
        
        [sg.Text('_'*30)],
        ]

    # ----- Full layout -----
    layout = [
        [sg.Column(first_col),
         sg.VSeperator(),
         sg.Column(seccond_col)]
    ]
    
    return(layout)


    


#%% Setup for run

path_to_data = '../../data/'
metadata_file = 'metadata.csv'
window = 10

meta = pd.read_csv(path_to_data+metadata_file)

# import each script as an ECM class item
data = []
for index,row in meta.iterrows():
    
    core = row['core']
    section = row['section']
    face = row['face']
    ACorDC = row['ACorDC']
    
    if core == 'alhic2302':
    #if core == 'pico2303':
        print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
        
        if section == '51_2':
            
            data_item = ECM(core,section,face,ACorDC)
            
            if max(data_item.depth) < 47 and min(data_item.depth)>3:
                data_item.smooth(window)
                data.append(data_item)

#launch gui
layout = make_gui()
window = sg.Window('Layerpicker GUI', layout, size=(1400, 900),
                   finalize=True, element_justification='center',
                   font='Helvetica 18')

#%% Run

qt = True
for d in data:
    
    print(" Running "+d.core+'-'+d.section+'-'+d.face+'-'+d.ACorDC)
    
    
    # read dataframe
    df = pd.read_csv(path_to_data+d.core+'/'+d.core+'-'+d.section+'-'+d.face+'-'+d.ACorDC+'.csv')
    if 'button_raw' in df.columns:
        print("Button_raw exists")
    else:
        df['button_raw'] = df['Button']
    
    # update current file
    window['-FILE-'].update(d.core+' '+d.section+' '+d.face+' '+d.ACorDC)

    # Loop through each track
    ycnt = 0
    for ycurr in d.y_vec:
        
        # update current track
        window['-TRACK-'].update(str(ycnt+1)+' of '+str(len(d.y_vec)))
    
            
        # set default plot limits
        ymin = min(d.depth_s) - (max(d.depth_s) -min(d.depth_s)) * 0.05
        d_window = 1.1 * (max(d.depth_s) -min(d.depth_s))
        
        # empty coords
        coords = 0

        # layer min and max vectors
        lmin = []
        lmax = []
        
        # Set plot true/false (catches if lmin>lmax)
        plot = True
        rescale = False
        update = True
        plt_but = False
        nextfile = False
        save = True
        
        # get initial x-axis bounds
        idx = np.logical_and(np.array(d.depth_s >= ymin),np.array(d.depth_s <= ymin+d_window))
        xmin = min(d.meas_s[idx])*0.95
        xmax = max(d.meas_s[idx])*1.05
        
        # add the plot to the window
        fig = makeplot(ymin,lmin,lmax,xmin,xmax,d_window,d,ycurr,df)
        tkcanvas = draw_figure(window['-CANVAS-'].TKCanvas, fig)
        
        while True:
            # read if button is pressed
            event, values = window.read(timeout=15)
        
            
            # activate button click
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            

            # check for quit
            if event in (sg.WIN_CLOSED, 'Quit'):
                qt = False
                break
            
            # check for move on to next file
            if event == '-NEXTTRACK-':
                break
            
            if event =='-NEXTFILE-':
                nextfile=True
                break
            
            if event =='-NEXTFILE_NOSAVE-':
                nextfile=True
                save=False
                break
            
            # check for down button
            elif event == '-DOWN-':
                ymin -= 0.02
                update = True
                
            # check for up button
            elif event == '-UP-':
                ymin += 0.02
                update = True
            
            # check for zoom in button
            elif event == '-+Z-':
                d_window *= (4/5)
                update = True
                
            # check for zoom out button
            elif event == '--Z-':
                d_window *= (5/4)
                update = True
            
            # if event is rescale
            elif event == '-RESCALE-':
                
                idx = np.logical_and(np.array(d.depth_s >= ymin),np.array(d.depth_s <= ymin+d_window))
                
                xmin = min(d.meas_s[idx])*0.95
                xmax = max(d.meas_s[idx])*1.05
                update = True
                #rescale = True
                
            # if event is rescale
            elif event == '-BUTTON-':
                if plt_but:
                    plt_but = False
                    window['-BUTTON-'].update('Show Button')
                else:
                    plt_but = True
                    window['-BUTTON-'].update('Hide Button')
                update = True
                
            
            # check for delete last
            elif event == '-DEL-':
                                
                # if there is a partial entry, only delete partial
                if len(lmax) == len(lmin):
                    lmax = lmax[:-1]
                
                lmin = lmin[:-1]
                
                # reset coords
                coords = 0
                
                update = True
                
            # check for new line clicked
            elif coords != 0:
                
                # if there is a new button click, then save it
                if (len(lmin)==0 or coords != lmin[-1] or plot == False) and (len(lmax)==0 or coords!= lmax[-1]):
        
                    # if we're waiting for a lmin
                    if len(lmax)==len(lmin):
                        lmin.append(coords)

                        window['-STATUS-'].update('Waiting for next Lmax')
                        
                    # if we're waiting for an lmax
                    else:
                        
                        # catch if the lmax is less than lmin, don't plot, loop
                        if coords > lmin[-1]:
                            lmax.append(coords)
                            plot = True
                            window['-STATUS-'].update('Waiting for next Lmin')
                        else:
                            plot = False
                            window['-STATUS-'].update('ERROR: Lmin>Lmax')
                            
                    if plot:
                        update = True
                        
            if update:
                delete_figure(tkcanvas)
                fig = makeplot(ymin,lmin,lmax,xmin,xmax,d_window,d,ycurr,df)
                tkcanvas = draw_figure(window['-CANVAS-'].TKCanvas, fig)
                window['-TBL-'].update(values=zip(lmin,lmax))
            
            # set update to false
            update = False
            rescale = False
        
        # update button column in dataframe
        if True:
            # zero values from htis track
            track_idx = df.index[df['Y_dimension(mm)'] == ycurr]
            df.loc[track_idx, 'Button'] = 0
            #loop through each lmin and lmax
            for i in range(min([len(lmax),len(lmin)])):
                idx = df.index[(df['True_depth(m)'] >= lmin[i]) & (df['True_depth(m)'] <= lmax[i]) & (df['Y_dimension(mm)'] == ycurr)] #& (df.index[df['Y_dimension(mm)'] == ycurr])]
                df.loc[idx, 'Button'] = 1
        
        # delete figure
        delete_figure(tkcanvas)
        
        # break out of loop if user requests next file
        if nextfile:
            break
        
        if not qt:
            break
    
        ycnt+=1
    
    
    # Adjust column
    if update_file and save:
        
        # save csv
        df.to_csv(path_to_data+d.core+'/'+d.core+'-'+d.section+'-'+d.face+'-'+d.ACorDC+'.csv',index=False)
    
    if not qt:
        break
    
window.close()
    
    
    