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

#origfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/S1/originals/'
#cropfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/S1/cropped/'
#pngfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/S1/png/'

year = '2018'

root = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/from_requests/max30daysep/nc/' + year + '/'
origfolder = root + 'orig/'
cropfolder = root + 'cropped/'
pngfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_TWITandTEIS/from_requests/max30daysep/png/' + year + '/cropped/'

# get all the nc files
origfiles = os.listdir(origfolder)

# get all the S1 img files to use as base images
imgfolder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'
imgfiles = [fn for fn in os.listdir(imgfolder) if fn.startswith('s1a')]
dates = [file[14:22] for file in imgfiles] #pull the date of the image from the filename

# turn crop toggle on if path to folder doesn't already exist
if os.path.exists(cropfolder):
    crop_nc = 0
else: 
    crop_nc = 1
    os.makedirs(cropfolder)
    os.makedirs(pngfolder)
    

#%% Crop all images

if crop_nc == 1:
    
    # Note: projection is in EPSG:3031 - WGS 84 / Antarctic Polar Stereographic
    # coordinates of bounding box. currently set to the extent of S1 images from Christian
    min_x = -1630295.6507051021326333
    min_y = -494928.3414874689187855 
    max_x = -1512646.2650982015766203
    max_y = -400908.0314721059403382
    
    for fn in origfiles:
        
        # only want the netCDF files
        if fn[-3:] == '.nc':
            
            try:
                nc = xr.open_dataset(origfolder + fn)
                
                # select only the subset that overlaps with bounding box above
                nc2 = nc.sel(y=slice(max_y,min_y), x=slice(min_x,max_x))
                
                # save the cropped netcdf!
                try:
                    nc2.to_netcdf(path=cropfolder + fn[:-3] +'_cropped' + fn[-3:], mode= "w")
                except:
                    print("not able to save file, might not overlap crop bounding box") 
                
            except:
                print("not able to open netCDF, might be empty") # 2024/03/12: urls to landsat velocity pairs are currently not working!
            
            

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
        
        ds_date = ds.img_pair_info.date_center[0:8]
        
        # find S1 base image with closest date
        # Calculate the absolute difference between each date and the given date
        differences = [abs(int(date) - int(ds_date)) for date in dates]

        # Find the index of the minimum difference
        ind = differences.index(min(differences))
        
        # get the image the corresponds to the date closest to the S1 image
        img = imgfolder + imgfiles[ind]
        
        # only if there's data
        if np.size(lat_) != 0 and np.size(lon_) != 0:
            # convert lat and lon to meshgrids
            lon, lat = np.meshgrid(lon_, lat_)
            
            # plotting parameters
            vmin = 10 #m/yr
            vmax = 5000 #m/yr
            title = ds_date
            n = 30
            
            savefig = 1
            outname = pngfolder + fn[:-3] + '.png'
            
            nq.nice_quiver(lon, lat, vx, vy, v, 
                           vmin, vmax, units = 'm/yr', n = n, img = img, title = title, 
                           savefig = savefig, outname = outname, 
                           cmap = 'magma_r')
            
            plt.close()
            