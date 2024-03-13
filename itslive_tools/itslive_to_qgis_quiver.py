#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 13:16:41 2024

@author: msharp
"""

import xarray as xr
import os

# %% Set up the paths

test_case = 0

if test_case == 1:  # an example file to test things out!

    inpath = "/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/"
    infn = "test_arrows.nc"
    
    goodfiles = infn
    
    outpath = "/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs_for_QGIS/"
    
    
elif test_case ==0:  # loop over all the data!
    
    #inpath = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs_shearmargin_90p_30days/'
    #outpath = "/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs_for_QGIS/velocity_pairs_shearmargin_90p_30days/"
    
    inpath = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs_for_QGIS/velocity_pairs/not_cropped/'
    outpath = "/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs_for_QGIS/velocity_pairs/"
    
    
    # get a list of all the files in the directory
    allfiles = os.listdir(inpath)
    
    # only want the netCDF files
    goodfiles = [file for file in allfiles if file.endswith('.nc') and not file.startswith('._') and not file.endswith('_cropped.nc')]

# %% Open the netcdf

for infn in goodfiles:
    
    # set the output filename to append "_cropped" but keep the extension ".nc"
    outfn = infn[:-3] + "_cropped" + ".nc"
   
    
    dc = xr.open_dataset(inpath+infn)
    
    # %% Crop the netcdf to bounding box
    
    # Note: projection is in EPSG:3031 - WGS 84 / Antarctic Polar Stereographic
    # coordinates of bounding box. currently set to the extent of S1 images from Christian
    min_x = -1630295.6507051021326333
    min_y = -494928.3414874689187855 
    max_x = -1512646.2650982015766203
    max_y = -400908.0314721059403382
    
    dc_cropped = dc.sel(y=slice(max_y,min_y), x=slice(min_x,max_x))
    
    # %% Rename velocity components for quiver plotting in QGIS
    
    dc2 = dc_cropped.rename(name_dict={'vx':'u-velocity component'})
    dc2 = dc2.rename(name_dict={'vy':'v-velocity component'})
    
    # %% save!
    
    try:
        dc2.to_netcdf(path=outpath+outfn, mode= "w")
    except:
        print("not able to save file. probably doesn't cover bounding box")
    
    if infn == goodfiles[len(goodfiles)-1]:
        
        print("Loop finished on file number: " + str(len(goodfiles)) + ", filename: " + outfn)
        