#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:07:38 2024

@author: msharp
"""

#### See Joanna Millstein's jupyter notebook for example!
#https://github.com/jdmillstein/strain_rates/blob/main/simple_cross_plots_crary.ipynb

import numpy as np
import matplotlib.pyplot as plt
import os, sys
from osgeo import gdal
import math as m
import matplotlib.colors as colors
import scipy.ndimage as nd
import xarray as xr
import matplotlib

from strain_codes import *

# %% download data

# downloaded example data of Crary ice rise into this folder 
#!wget https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0754.001/1996.01.01/antarctic_ice_vel_phase_map_v01.nc 

# %% open velocity dataset and split into vx and vy components

ds = xr.open_dataset(('antarctic_ice_vel_phase_map_v01.nc'))
vx_ = np.asarray(ds.VX)
vy_ = np.asarray(ds.VY)
lat_ = np.asarray(ds.lat)
lon_ = np.asarray(ds.lon)

# %% crop to region and visualize

#use plt.imshow(vx) to crop dimensions to some region that is both computationally efficient and small enough!

vx = vx_[7300:8100,5000:6600]
vy = vy_[7300:8100,5000:6600]
lat = np.flip((lat_[7300:8100,5000:6600]).T)
lon = np.flip((lon_[7300:8100,5000:6600]).T)-360

plt.imshow(vx, vmin = 10, vmax = 300, cmap = 'magma',
           extent = [lon[-1,0], lon[0,0],  lat[0,0], lat[0,-1]], 
           aspect = 10)

# %% Compute strain rates

#compute strain rates *contained in the dictionary crary* and the strain rate tensor *DR* from Bryan Riel's iceutils
#(this code is a slightly simplified version to mess with the window size)
# strain rates spit out in values of per year

crary, DR = compute_strains(vx, vy, 450, -450, rotate=True, grad_method='numpy',inpaint=False)

#compute the eigenvalues and eigenvectors
eigvalsR, eigvecsR = np.linalg.eigh(DR.T)

princ1stR = eigvalsR[:,1].reshape(vx.shape)
princ2ndR = eigvalsR[:,0].reshape(vx.shape)
ev01R = eigvecsR[:,0,1].reshape(vx.shape)
ev00R = eigvecsR[:,0,0].reshape(vx.shape)
ev11R = eigvecsR[:,1,1].reshape(vx.shape)
ev10R = eigvecsR[:,1,0].reshape(vx.shape)
eigvals, eigvecs = np.linalg.eigh(DR.T)

# %% Build colourmap

# colormap where > 0 is blue (extension) and < 0 is red (compression), can remove 

bwr = matplotlib.colormaps['bwr']
newcolors = bwr(np.linspace(0, 1, 256))

bone =  matplotlib.colormaps['bone']
bone = bone(np.linspace(0, 1, 100))

newcolors2 = newcolors[0:100, :] 
greys = bone[70:,:]; 
newcolors3 = newcolors[156:, :]
newcmp = colors.ListedColormap(np.vstack((newcolors3[::-1], np.tile(greys[-1],[10,1]), newcolors2[::-1])))
newcmp

# %% Plot principal components!

# plot cropped region around Brunt's north rift in lat/lon S. polar stereo   [ lon[0,0],  lat[0,0], lat[0,-1]]
lon = np.linspace(-1*lon[-1,0], -1*lon[0,0], len(vx[1,:]))
lon = -1*(lon)
lat = np.linspace(lat[0,0],lat[0,-1],len(vx[:,1]))


la,lo= np.meshgrid(lon,lat)

fig, axs = plt.subplots(ncols=1,dpi=205)

i = 60
j = 60

velocity = np.sqrt(vx**2 + vx**2)

sqrt = ((np.sqrt(vx[::i,::j]**2 + vy[::i,::j]**2)))

# velocity in meters per year
axs.imshow(nd.median_filter(velocity,5), vmin = 0, vmax = 3e2, cmap='magma', alpha=0.25,
           extent = [lon[0], lon[-1], lat[-1], lat[0]], aspect=10)

# first principal component 
axs.quiver(la[::i,::j],lo[::i,::j],
                ev00R[::i,::j], 
                ev01R[::i,::j], np.ma.masked_where(sqrt > 40000.,princ1stR[::i,::j]),
                pivot='middle',headaxislength=0, headlength=0,angles='uv',
                cmap=newcmp,norm=colors.TwoSlopeNorm(vmin=-.01, vcenter=0,vmax=.01))

# second principal component of strain
axs.quiver(la[::i,::j],lo[::i,::j],
                ev10R[::i,::j],
                ev11R[::i,::j], np.ma.masked_where(sqrt > 40000.,princ2ndR[::i,::j]),
                pivot='middle',headaxislength=0, headlength=0,angles='uv',
                cmap=newcmp,norm=colors.TwoSlopeNorm(vmin=-.01, vcenter=0,vmax=.01))

#axs.axis('off')
plt.suptitle("Principal Components of Strain")
plt.title('brighter color indicates larger eigenvalue (P1 or P2). Red is compression, Blue extension')
plt.show()

# %% Plot effective, longitudinal, transverse, and shear strain rate components

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

# Plotting each dataset in the grid
im = axs[0, 0].imshow(crary['effective'], cmap='PuBu', vmin=0, vmax = 0.1)
cbar = fig.colorbar(im, ax=axs[0, 0]); cbar.set_label('$yr^{-1}$')
axs[0, 0].set_title('effective strain rate')

im = axs[0, 1].imshow(crary['e_xx'], cmap='bwr_r', vmin = -.05, vmax = .05)
cbar = fig.colorbar(im, ax=axs[0, 1]); cbar.set_label('$yr^{-1}$')
axs[0, 1].set_title('$\epsilon_{xx} \;$ (along-flow)')

im = axs[1, 0].imshow(crary['e_yy'], cmap='bwr_r', vmin = -.05, vmax = .05)
cbar = fig.colorbar(im, ax=axs[1, 0]); cbar.set_label('$yr^{-1}$')
axs[1, 0].set_title('$\epsilon_{yy}\;$ (across-flow)')

im = axs[1, 1].imshow(crary['e_xy'], cmap='PuBu', vmin = 0, vmax = 0.1)
cbar = fig.colorbar(im, ax=axs[1, 1]); cbar.set_label('$yr^{-1}$')
axs[1, 1].set_title('shear strain rate')


plt.tight_layout()
plt.show()