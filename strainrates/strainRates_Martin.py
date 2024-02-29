#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 10:13:05 2023

@author: msharp

"""

### This code is adapted from Martin Truffer!

import xarray as xr
import numpy as np
import matplotlib as plt


# %% Load dataset
#recon = xr.open_dataset('/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs/S1A_IW_SLC__1SSH_20150701T043600_20150701T043627_006612_008D15_59D2_X_S1A_IW_SLC__1SSH_20150713T043548_20150713T043616_006787_0091F9_6500_G0120V02_P062.nc')
recon = xr.open_dataset('/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs_for_QGIS/velocity_pairs/S1A_IW_SLC__1SSH_20150701T043600_20150701T043627_006612_008D15_59D2_X_S1A_IW_SLC__1SSH_20150713T043548_20150713T043616_006787_0091F9_6500_G0120V02_P062_cropped.nc')

vx = recon.vx
vy = recon.vy
v = recon.v
y = recon.y
x = recon.x

# %% convert to numpy and drop nans

v = v.to_numpy()

vx = vx.to_numpy()
vx = vx[row_mask, col_mask]

vy = vy.to_numpy()
vy = vy[~np.isnan(v)]

y = y.to_numpy()
y = y[~np.isnan(v)]

x = x.to_numpy()
x = x[~np.isnan(v)]

v = v[~np.isnan(v)]


# %% Calculate the 2D mean of each time slice falling in the same year

red = 20

# slice from start to end in "red" number of chunks
vx_red = vx[::red,::red].to_numpy()
vy_red = vy[::red,::red].to_numpy()
x_red = x[::red].to_numpy()
y_red = y[::red].to_numpy()

dx = np.max(np.diff(x_red))
dy = np.max(np.diff(y_red))

# %% calculate derivatives and strain rates

dvxdx = (vx_red[1:-1,2:]-vx_red[1:-1,:-2])/(2*dx)
dvxdy = (vx_red[2:,1:-1]-vx_red[:-2,1:-1])/(2*dy)
dvydx = (vy_red[1:-1,2:]-vy_red[1:-1,:-2])/(2*dx)
dvydy = (vy_red[2:,1:-1]-vy_red[:-2,1:-1])/(2*dy)

exy = 0.5*(dvxdy+dvydx)

e1 = 0.5*(dvxdx+dvydy) + 0.5*np.sqrt((dvxdx-dvydy)**2+4*exy**2)
e2 = 0.5*(dvxdx+dvydy) - 0.5*np.sqrt((dvxdx-dvydy)**2+4*exy**2)

# %% principal directions

n1 = np.sqrt(exy**2+(dvxdx-e1)**2)
pd1x = exy/n1
pd1y = (dvxdx-e1)/n1

n2 = np.sqrt(exy**2+(dvxdx-e2)**2)
pd2x = exy/n2
pd2y = (dvxdx-e2)/n2

def plot_strainrates(x_red, y_red, e1,e2,pd1x,pd1y,pd2x,pd2y):
    scale=10000
    i=0
    for xi in x_red[1:-1]:
        j=0
        for yi in y_red[1:-1]:
            if e1[j,i]>0:
                clr='r'
            else:
                clr='b'
            plt.plot([xi-pd1x[j,i]*scale*e1[j,i], xi+pd1x[j,i]*scale*e1[j,i]],
                     [yi-pd1y[j,i]*scale*e1[j,i], yi+pd1y[j,i]*scale*e1[j,i]], color=clr)
            if e2[j,i]>0:
                clr='r'
            else:
                clr='b'
            plt.plot([xi-pd2x[j,i]*scale*e2[j,i], xi+pd2x[j,i]*scale*e2[j,i]],
                     [yi-pd2y[j,i]*scale*e2[j,i], yi+pd2y[j,i]*scale*e2[j,i]], color=clr)
            j+=1
        i+=1
    plt.axis('equal')
       

plt.figure(figsize=(16,12))
plt.pcolormesh(x,y,v, shading='nearest', vmax=3000, alpha=1, cmap='cool')
clb=plt.colorbar()
plot_strainrates(x_red, y_red, e1,e2,pd1x,pd1y,pd2x,pd2y)