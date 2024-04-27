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

# gui
import PySimpleGUI as sg

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM

global fig

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
    
#%% return coorinates when I click on the graph
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print (f'x = {ix}, y = {iy}')

    global coords
    coords = iy
    
    fig.canvas.mpl_disconnect(cid)

    return iy


# make plot
def makeplot(ymin,lmin,lmax,xmin,xmax,d_window,data,currtrack):


    yspc = data.y_vec[1]-data.y_vec[0]
    cmap = matplotlib.colormaps.get_cmap('coolwarm')
    
    fig, ax = plt.subplots(1, 2,figsize=(10, 7), dpi=100)
    
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
        ind2 = ind * data.button_s==1
        ax[1].plot(data.meas_s[ind2],data.depth_s[ind2],'k.',markersize=linewidth)
    
    cnt = 0
    for x in lmin:

        ax[1].plot([min(data.meas_s), max(data.meas_s)],[x,x],'k')
        ax[0].plot([currtrack-yspc/2, currtrack+yspc/2],[x,x],'k')

        cnt +=1
        
    cnt = 0
    for x in lmax:
        
        ax[1].plot([min(data.meas_s), max(data.meas_s)],[x,x],'g')
        ax[1].add_patch(Rectangle((min(data.meas_s),lmin[cnt]),max(data.meas_s)-min(data.meas_s),x-lmin[cnt],facecolor=(1, 0, 0, 0.2)))

        ax[0].plot([currtrack-yspc/2, currtrack+yspc/2],[x,x],'g')
        ax[0].add_patch(Rectangle((currtrack-yspc/2,lmin[cnt]),yspc,x-lmin[cnt],facecolor=(1, 0, 0, 0.2)))


        cnt +=1

    ax[1].set_xlim(xmin,xmax)
    ax[1].set_ylim(ymin,ymin+d_window)
    ax[1].invert_yaxis()
    ax[1].set_ylabel('Depth (m)')
    ax[1].set_xlabel('Conductivity (amps)')
    ax[1].set_title('Conductivity Plots')
    
    ax[0].set_ylim(ymin+d_window,ymin)
    ax[0].invert_xaxis()
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Conductivity (amps)')
    ax[0].set_title('Conductivity Plots')
    
    fig.tight_layout()
    
    return(fig)

#%% make GUI

def make_gui():
    
    # Create GUI
    first_col = [
       [sg.Text('Plot test')],
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
        [sg.Text('Current File'),sg.Text(size=(25, 1), key='-FILE-', pad=(1, 1))],
        [sg.Button(button_text='Go To Next Track',key = '-NEXTTRACK-')],
        [sg.Button(button_text='Go To Next File (save current work)',key = '-NEXTFILE-')],
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


    


#%% Test

path_to_data = '../../data/'
metadata = 'metadata.csv'

    
data = [ECM('alhic2302','38','t','AC'),ECM('alhic2302','39','t','AC')]

for d in data:
    d.smooth(10)
    
    d.norm_all()

#launch gui
layout = make_gui()
window = sg.Window('Layerpicker GUI', layout, size=(1400, 900),
                   finalize=True, element_justification='center',
                   font='Helvetica 18')

qt = True
for d in data:
    
    for ycurr in d.y_vec:
    
        if qt:
            
            # set default plot limits
            ymin = min(d.depth_s)
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
            
            # get initial x-axis bounds
            idx = np.logical_and(np.array(d.depth_s >= ymin),np.array(d.depth_s <= ymin+d_window))
            xmin = min(d.meas_s[idx])*0.95
            xmax = max(d.meas_s[idx])*1.05
            
            # add the plot to the window
            fig = makeplot(ymin,lmin,lmax,xmin,xmax,d_window,d,ycurr)
            tkcanvas = draw_figure(window['-CANVAS-'].TKCanvas, fig)
            
            # update current file
            window['-FILE-'].update(d.core+' '+d.section+' '+d.face+' '+d.ACorDC)
            
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
                    
                    print('Deleting')
                    
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
                            print('setting update to True')
                            
                if update:
                    delete_figure(tkcanvas)
                    fig = makeplot(ymin,lmin,lmax,xmin,xmax,d_window,d,ycurr)
                    tkcanvas = draw_figure(window['-CANVAS-'].TKCanvas, fig)
                    window['-TBL-'].update(values=zip(lmin,lmax))
                
                # set update to false
                update = False
                rescale = False
            
        else:
            break
        
        delete_figure(tkcanvas)
        if nextfile:
            break
        
    
window.close()
    
    
    