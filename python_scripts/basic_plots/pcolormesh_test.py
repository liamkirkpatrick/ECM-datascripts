#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 11:55:20 2024

@author: Liam
"""

import numpy as np
import matplotlib.pyplot as plt

# my functions/classes
import sys
sys.path.append("../core_scripts/")
from ECMclass import ECM


path_to_data = '../../data/'
metadata = 'metadata.csv'

data = ECM('alhic2302','38','t','AC')

data.smooth(10)

# Example vectors
x = data.y_s #-(data.y_vec[1]-data.y_vec[0])/2  # x dimension points
y = data.depth_s    # y dimension points
z = data.meas_s # Value at each point


# Find the unique grid points
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
plt.pcolormesh(X, Y, Z, shading='auto')
plt.colorbar()  # Show color scale
plt.show()
