#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 17:26:41 2024

@author: msharp
"""

import nice_quiver as nq
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# %% Paths

origfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/S1/originals/'
cropfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/S1/cropped/'
pngfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/S1/png/'

origfiles = os.listdir(origfolder)

img = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/s1a-iw-grd-hh-20210530t043630-20210530t043655-038112-047f7a-001_ps_ROI.tiff'

# %% Crop all images. ONLY uncomment if you want to recrop and save!

# # Note: projection is in EPSG:3031 - WGS 84 / Antarctic Polar Stereographic
# # coordinates of bounding box. currently set to the extent of S1 images from Christian
# min_x = -1630295.6507051021326333
# min_y = -494928.3414874689187855 
# max_x = -1512646.2650982015766203
# max_y = -400908.0314721059403382

# for fn in origfiles:
    
#     # only want the netCDF files
#     if fn[-3:] == '.nc':
#         nc = xr.open_dataset(origfolder + fn)
        
#         # select only the subset that overlaps with bounding box above
#         nc2 = nc.sel(y=slice(max_y,min_y), x=slice(min_x,max_x))
        
#         # save the cropped netcdf!
#         try:
#             nc2.to_netcdf(path=cropfolder + fn[:-3] +'_cropped' + fn[-3:], mode= "w")
#         except:
#             print("not able to save file. probably doesn't cover bounding box")
 

# %% Plots

cropfiles = os.listdir(cropfolder)

for fn in cropfiles:
    
    # only want the netCDF files
    if fn[-3:] == '.nc':
        
        # open the netCDF dataset
        ds = xr.open_dataset(cropfolder + fn)
        
        # create numpy arrays
        vx = np.asarray(ds.vx)
        vy = np.asarray(ds.vy)
        v = np.asarray(ds.v)
        lat_ = np.asarray(ds.y)
        lon_ = np.asarray(ds.x)
        
        # only if there's data
        if np.size(lat_) != 0 and np.size(lon_) != 0:
            # convert lat and lon to meshgrids
            lon, lat = np.meshgrid(lon_, lat_)
            
            # plotting parameters
            vmin = 10 #m/yr
            vmax = 5000 #m/yr
            title = fn[17:25]
            n = 30
            
            savefig = 1
            outname = pngfolder + fn[:-3] + '.png'
            
            nq.nice_quiver(lon, lat, vx, vy, v, 
                           vmin, vmax, units = 'm/yr', n = n, img = img, title = title, 
                           savefig = savefig, outname = outname, 
                           cmap = 'magma_r')
            
            plt.close()
            