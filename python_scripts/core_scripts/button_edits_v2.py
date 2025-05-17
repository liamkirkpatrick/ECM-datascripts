import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg

from ECMclass import ECM

class LayerPicker:
    def __init__(self, data_item, window):
        self.d = data_item
        self.df = pd.read_csv(f"{path_to_data}{self.d.core}/{self.d.core}-{self.d.section}-{self.d.face}-{self.d.ACorDC}.csv")
        if 'button_raw' not in self.df:
            self.df['button_raw'] = self.df['Button']
        self.lmin, self.lmax = [], []
        self.curr_track_idx = 0
        self.ymin = self.d.depth_s.min() - 0.05 * (self.d.depth_s.ptp())
        self.d_window = 1.1 * self.d.depth_s.ptp()
        self._build_figure()
        self._draw_initial()
        self._connect_events()
        self.window = window
        self._blit_background = None

    def _build_figure(self):
        self.fig, self.axes = plt.subplots(1, 2, figsize=(10,7), dpi=100)
        # precompute meshgrid once
        x_unique = np.unique(self.d.y_s)
        y_unique = np.unique(self.d.depth_s)
        X, Y = np.meshgrid(x_unique, y_unique)
        Z = np.full_like(X, np.nan, dtype=float)
        for xi, yi, zi in zip(self.d.y_s, self.d.depth_s, self.d.meas_s):
            ix = np.searchsorted(x_unique, xi)
            iy = np.searchsorted(y_unique, yi)
            Z[iy, ix] = zi
        self.mesh = self.axes[0].pcolormesh(X, Y, Z, shading='auto', cmap=matplotlib.colormaps['Spectral'])
        # lines for each track
        self.lines = []
        cmap = matplotlib.colormaps['coolwarm']
        for i, yval in enumerate(self.d.y_vec):
            line, = self.axes[1].plot([], [], lw=2, color=cmap(i/len(self.d.y_vec)))
            self.lines.append(line)
        # artist for button points
        self.btn_scatter = self.axes[1].scatter([], [], s=20, c='k')
        # patches for lmin/lmax rectangles
        self.rects = []
        # highlight current track
        self.track_rect = Rectangle((0,0), 1,1, edgecolor='black', facecolor='none')
        self.axes[0].add_patch(self.track_rect)

        for ax in self.axes:
            ax.set_xlabel('Conductivity (amps)')
            ax.set_ylabel('Depth (m)')
        self.axes[0].set_title('Top View')
        self.axes[1].set_title('Individual Curves')
        self.fig.tight_layout()

    def _draw_initial(self):
        self._update_plot()
        # blit background
        canvas = self.canvas_widget
        self.canvas_widget.draw()
        self._blit_background = canvas.copy_from_bbox(self.fig.bbox)

    def _connect_events(self):
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)

    def _on_click(self, event):
        if event.inaxes != self.axes[1]: return
        y_click = event.ydata
        if len(self.lmin) == len(self.lmax):
            self.lmin.append(y_click)
            self.window['-STATUS-'].update('Waiting for next Lmax')
        else:
            if y_click > self.lmin[-1]:
                self.lmax.append(y_click)
                self.window['-STATUS-'].update('Waiting for next Lmin')
            else:
                self.window['-STATUS-'].update('ERROR: Lmin>Lmax')
        self._update_plot(blit=True)

    def _update_plot(self, blit=False):
        # update current track label
        curr_y = self.d.y_vec[self.curr_track_idx]
        self.window['-TRACK-'].update(f"{self.curr_track_idx+1} of {len(self.d.y_vec)}")
        # update line data
        for line, yval in zip(self.lines, self.d.y_vec):
            mask = self.d.y_s == yval
            line.set_data(self.d.meas_s[mask], self.d.depth_s[mask])
            line.set_linewidth(7 if yval==curr_y else 2)
        # update button scatter
        mask_btn = (self.df['Y_dimension(mm)']==curr_y) & (self.df['Button']==1)
        xs = self.df.loc[mask_btn, 'Conductivity']
        ys = self.df.loc[mask_btn, 'True_depth(m)']
        self.btn_scatter.set_offsets(np.column_stack([xs, ys]))
        # update mesh and view limits
        self.axes[1].set_ylim(self.ymin, self.ymin + self.d_window)
        idx = (self.d.depth_s>=self.ymin)&(self.d.depth_s<=self.ymin+self.d_window)
        xm, xM = self.d.meas_s[idx].min()*0.95, self.d.meas_s[idx].max()*1.05
        self.axes[1].set_xlim(xm, xM)
        self.axes[0].set_ylim(self.ymin, self.ymin + self.d_window)
        self.axes[0].set_xlim(curr_y - self.d.y_vec.ptp()/len(self.d.y_vec),
                              curr_y + self.d.y_vec.ptp()/len(self.d.y_vec))
        # redraw rectangles
        for rect in self.rects:
            rect.remove()
        self.rects.clear()
        for ymin, ymax in zip(self.lmin, self.lmax):
            rect = Rectangle((xm, ymin), xM-xm, ymax-ymin, facecolor=(1,0,0,0.2))
            self.axes[1].add_patch(rect)
            self.rects.append(rect)
        # track rect
        self.track_rect.set_xy((curr_y - self.d.y_vec.ptp()/len(self.d.y_vec)/2, self.ymin))
        self.track_rect.set_height(self.d_window)
        self.track_rect.set_width(self.d.y_vec.ptp()/len(self.d.y_vec))

        if blit:
            canvas = self.canvas_widget
            canvas.restore_region(self._blit_background)
            for ax in self.axes:
                ax.draw_artist(self.mesh)
                for line in self.lines:
                    ax.draw_artist(line)
                for rect in self.rects + [self.track_rect]:
                    ax.draw_artist(rect)
            canvas.blit(self.fig.bbox)
        else:
            self.canvas_widget.draw()

    def next_track(self):
        # save current track button settings back to df
        curr_y = self.d.y_vec[self.curr_track_idx]
        self.df.loc[self.df['Y_dimension(mm)']==curr_y, 'Button'] = 0
        for ymin, ymax in zip(self.lmin, self.lmax):
            mask = ((self.df['True_depth(m)']>=ymin)&
                    (self.df['True_depth(m)']<=ymax)&
                    (self.df['Y_dimension(mm)']==curr_y))
            self.df.loc[mask, 'Button'] = 1

        # reset for new track
        self.curr_track_idx += 1
        self.lmin.clear()
        self.lmax.clear()
        if self.curr_track_idx >= len(self.d.y_vec):
            return False
        self._update_plot()
        return True

    def save_and_close(self):
        self.df.to_csv(f"{path_to_data}{self.d.core}/{self.d.core}-{self.d.section}-{self.d.face}-{self.d.ACorDC}.csv", index=False)

    @property
    def canvas_widget(self):
        return self._canvas

    @canvas_widget.setter
    def canvas_widget(self, tkagg):
        self._canvas = tkagg

# --- main GUI loop ---

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
    to_run = ['158']

    if core == 'alhic2416':
        
        print("Reading "+core+", section "+section+'-'+face+'-'+ACorDC)
        
        data_item = ECM(core,section,face,ACorDC)
        
        #if max(data_item.depth) < 47 and min(data_item.depth)>3:
        data_item.smooth(window)
        data.append(data_item)

layout = make_gui()  # same as your layout fn
win = sg.Window('Layerpicker', layout, finalize=True, size=(1400,900))
canvas_elem = win['-CANVAS-'].TKCanvas

for d in data_objs:
    picker = LayerPicker(d, win)
    picker.canvas_widget = FigureCanvasTkAgg(picker.fig, canvas_elem)
    picker.canvas_widget.get_tk_widget().pack(side='top', fill='both', expand=1)

    while True:
        event, _ = win.read(timeout=50)
        if event in (sg.WIN_CLOSED, 'Quit'):
            raise SystemExit
        elif event == '-NEXTTRACK-':
            if not picker.next_track():
                picker.save_and_close()
                break
        # handle zoom/up/down/NEXTFILE etc. by adjusting picker.ymin or picker.d_window
        # then call picker._update_plot(blit=True)

win.close()
