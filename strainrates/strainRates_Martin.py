#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 10:13:05 2023

@author: msharp, adopted from Martin Truffer

"""

#%% This code is adapted from Martin Truffer!

import xarray as xr
import rioxarray as rioxr
import numpy as np
import rasterio

import matplotlib
from matplotlib import pyplot as plt
#matplotlib.use('Agg') # Non-interactive backend for file saving
matplotlib.use("inline")

import os
import geopandas as gpd

# %% toggles
itslive = 0
adrian = 1
load_img = 1
crop_vel = 0

compare_dates = 0

gash = 0
channel = 0
daggers = 0

specify_pairs = 0

principal_directions = 1
divergence = 0


savepng = 0

# %% helper functions

def plot_strainrates(x_red, y_red, e1,e2, pd1x,pd1y,pd2x,pd2y, col1, col2):
    
    print('ploting...')
    
    fig1, ax1 = plt.subplots(ncols=1,dpi=205)
    #ax1.set_aspect('equal')
    #fig2, ax2 = plt.subplots(ncols=1,dpi=205)
    
    # if col1 == 'xkcd:blurple':
    #     ax1.imshow(img, vmin= 0, vmax = 400, cmap = 'gray', extent = extent) 
    
    ax1.imshow(img, vmin= 0, vmax = 400, cmap = 'gray', extent = extent) 
    #ax1.set_xlim(-1.60494e6,-1.56003e6)
    #ax1.set_ylim(-4.6592e5, -4.2995e5)
    #ax2.imshow(img, vmin= 0, vmax = 400, cmap = 'gray', extent = extent)
    
    
    #ax1.pcolormesh(x,y,v, shading='nearest', vmax=3000, alpha=0.3, cmap='magma_r')
    #ax2.pcolormesh(x,y,v, shading='nearest', vmax=3000, alpha=0.3, cmap='magma_r')
    
    # full extent
    # scale=2500
    # threshold=1
    
    # i think these settings are for Daggers?
    # scale=20000
    # threshold=0.1
    
    scale = 1000
    threshold = 50
    
    
    i=0
    
    for xi in x_red[1:-1]:
        j=0
        
        for yi in y_red[1:-1]:
            
            # Confused about signs???
            # Martin had negative = red (compression?)
            # but other papers have:
            # compression is negative (blue)
            # extension is positive (red)
            
           # plot smallest strain rate component fist
            if e2[j,i]<0:
                #clr='b'
                clr = col2
                       
            else:
                #clr='r'
                clr = col1
                
            if abs(e2[j, i]) < threshold and abs(e1[j, i]) < threshold:
                ax1.plot([xi-pd2x[j,i]*scale*e2[j,i], xi+pd2x[j,i]*scale*e2[j,i]],
                         [yi-pd2y[j,i]*scale*e2[j,i], yi+pd2y[j,i]*scale*e2[j,i]], color=clr, linewidth = 0.8)


            # plot largest strain rate component on top
            if e1[j,i]<0: 
                #clr='r' # Martin had the colours this way
                #clr = 'b'
                clr = col2
                
            else:
                #clr='r'
                # clr = col1
                clr = col1
            
            if abs(e2[j, i]) < threshold and abs(e1[j, i]) < threshold:
                ax1.plot([xi-pd1x[j,i]*scale*e1[j,i], xi+pd1x[j,i]*scale*e1[j,i]],
                         [yi-pd1y[j,i]*scale*e1[j,i], yi+pd1y[j,i]*scale*e1[j,i]], color=clr, linewidth = 0.8)
            
                
               
            
            j+=1
        i+=1
    
    return fig1, ax1
    
    

def find_closest_image(date, folder):
    
    # make sure the input date is an integer (yyyymmdd)
    date = int(date)
    
    img_list = [im for im in os.listdir(folder) if im.endswith('_ROI.tiff') and not im.startswith('._')]
    
    date_diffs = []
    for fn in img_list:
        img_date = int(fn[14:22])
        
        date_diffs.append(abs(img_date-date))
        
    closest_idx = np.argmin(date_diffs)
    closest_img = img_list[closest_idx]
    
    return closest_img
        
    

# %% set up paths

img_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'

img_file2017 = 's1a-iw-grd-hh-20171217t043611-20171217t043636-019737-0218ee-001_ps_ROI.tiff'
img_file2018 = 's1a-iw-grd-hh-20180510t043611-20180510t043636-021837-025b44-001_ps_ROI.tiff'
img_file2019 = 's1a-iw-grd-hh-20190529t043618-20190529t043643-027437-031864-001_ps_ROI.tiff'
img_file2021 = 's1a-iw-grd-hh-20210530t043630-20210530t043655-038112-047f7a-001_ps_ROI.tiff'
img_file2022 = 's1a-iw-grd-hh-20220113t043634-20220113t043659-041437-04ed4e-001_ps_ROI.tiff'


# %%
if itslive:
    
    # # path to ITS_LIVE velocity pairs
    folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs_for_QGIS/velocity_pairs_shearmargin_90p_30days/'
    
    # used the cropped netCDFs
    filename2017 = 'S1A_IW_SLC__1SSH_20171216T084647_20171216T084714_019725_02188D_D67E_X_S1B_IW_SLC__1SSH_20171222T084605_20171222T084632_008829_00FB9B_9FF0_G0120V02_P091_cropped.nc'
    filename2018 = 'S1B_IW_SLC__1SSH_20180503T084605_20180503T084632_010754_013A64_E1E9_X_S1A_IW_SLC__1SSH_20180509T084647_20180509T084714_021825_025AE3_8978_G0120V02_P091_cropped.nc'
    filename2019 = 'S1B_IW_SLC__1SSH_20190522T084612_20190522T084639_016354_01EC80_4F45_X_S1A_IW_SLC__1SSH_20190528T084654_20190528T084721_027425_031802_570D_G0120V02_P091_cropped.nc'
    filename2021 = 'S1B_IW_SLC__1SSH_20210524T043547_20210524T043615_027041_033B0D_856F_X_S1A_IW_SLC__1SSH_20210530T043629_20210530T043657_038112_047F7A_3562_G0120V02_P092_cropped.nc'
    filename2022 = 'S2B_MSIL1C_20220111T151259_N0301_R139_T13CDS_20220111T181207_X_S2B_MSIL1C_20220121T151259_N0301_R139_T13CDS_20220121T181438_G0120V02_P096_cropped.nc'
    
    
    # choose which dataset
    year = 2017
    
    if year == 2017:
        filename = filename2017
        img_file = img_file2017
        
    elif year == 2018:
        filename = filename2018
        img_file = img_file2018
        
    elif year == 2019:
        filename = filename2019
        img_file = img_file2019
            
    elif year == 2021:
        filename = filename2021
        img_file = img_file2021
        
    elif year == 2022:
        filename = filename2022
        img_file = img_file2022
    
    #  Load ITSLIVE dataset
    #recon = xr.open_dataset('/Volumes/Sandisk4TB/PhD_MS/TARSAN/ITS_LIVE/netCDFs/velocity_pairs/S1A_IW_SLC__1SSH_20150701T043600_20150701T043627_006612_008D15_59D2_X_S1A_IW_SLC__1SSH_20150713T043548_20150713T043616_006787_0091F9_6500_G0120V02_P062.nc')
    recon = xr.open_dataset(folder + filename)
    
    # undo the change in variable name that was done in itslive_to_qgis_quiver.py when netcdfs were cropped
    recon = recon.rename(name_dict={'u-velocity component':'vx'})
    recon = recon.rename(name_dict={'v-velocity component':'vy'})
    
    vx = recon.vx
    vy = recon.vy
    v = recon.v
    y = recon.y
    x = recon.x

# %%
if adrian:
    
    #fig1, ax1 = plt.subplots(ncols=1,dpi=205)
    
    # Use Adrian's velocity mosaics


    # only do 2017 vs 2022
    x2017_04 = 'mean_mosaic.201704_201704.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.x.tif'
    y2017_04 = 'mean_mosaic.201704_201704.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.y.tif'
    
    x2018_08 = 'mean_mosaic.201808_201808.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.x.tif'
    y2018_08 = 'mean_mosaic.201808_201808.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.y.tif'
    
    x2019_11 = 'mean_mosaic.201911_201911.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.x.tif'
    y2019_11 = 'mean_mosaic.201911_201911.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.y.tif'
    
    x2020_01 = 'mean_mosaic.202001_202001.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.x.tif'
    y2020_01 = 'mean_mosaic.202001_202001.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.y.tif'
    
    
    x2021_08 = 'mean_mosaic.202108_202108.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.x.tif'
    y2021_08 = 'mean_mosaic.202108_202108.6or12-day.1to3-pairs.tracks-007.mosaics-10.track_filt_detide.mean.y.tif'
    
    x2021_12 = 'mean_mosaic.202112_202112.6or12-day.1to3-pairs.tracks-007.mosaics-8.track_filt_detide.mean.x.tif'
    y2021_12 = 'mean_mosaic.202112_202112.6or12-day.1to3-pairs.tracks-007.mosaics-8.track_filt_detide.mean.y.tif'
    
    
    if gash == 1:
        
        vel_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Velocity_and_Strain/Velocity_Fields_from_Adrian/Velocity_Monthly_2015_2023_Rishi_Adrian/'
        
        x_list = [x2018_08, x2021_08]
        y_list = [y2018_08, y2021_08]
        
    elif channel == 1:
        
        vel_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Velocity_and_Strain/Velocity_Fields_from_Adrian/Velocity_Monthly_2015_2023_Rishi_Adrian/'
        
        x_list = [x2020_01, x2021_12]
        y_list = [y2020_01, y2021_12]
        
    elif specify_pairs == 1:
        
        x_list = [x2017_04, x2021_08]
        y_list = [y2017_04, y2021_08]
        
    elif daggers == 1:
        
        vel_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Velocity_and_Strain/Velocity_Fields_from_Adrian/Velocity_Monthly_2016_2022_Adrian/'

        # get all the x grids
        x_list = [fn for fn in os.listdir(vel_folder) if fn.endswith('track_detide.x.tif') and '6-day' in fn and not fn.startswith('._')]
        # # get all the y grids
        y_list = [fn for fn in os.listdir(vel_folder) if fn.endswith('track_detide.y.tif') and '6-day' in fn and not fn.startswith('._')]
        
    else:
        
        # vel_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Velocity_and_Strain/Velocity_Fields_from_Adrian/Velocity_Monthly_2016_2022_Adrian/'
        vel_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Velocity_and_Strain/Velocity_Fields_from_Adrian/TEIS2016to2025/'
        
        # get all the x grids
        x_list = [fn for fn in os.listdir(vel_folder) if fn.endswith('.x.tif') and not fn.startswith('._')]
        # # get all the y grids
        y_list = [fn for fn in os.listdir(vel_folder) if fn.endswith('.y.tif') and not fn.startswith('._')]
        
 # %%   

    for ind,vel in enumerate(x_list):
        
        date = int(vel[12:18] + '15')
        # date = int(vel[7:13] + '15')
        print('working on: ' + vel)
        
        img_file = find_closest_image(date, img_folder)
        
        if compare_dates == 1:
            
            if date == 20170415:
                load_img = 0
                col1 = 'xkcd:marigold'
                col2 = col1
            elif date == 20210815:
                load_img = 1
                col1 = 'xkcd:blurple'
                col2 = col1
               
    
        # load the velocity data
        vx = rioxr.open_rasterio(vel_folder + x_list[ind])
        vy = rioxr.open_rasterio(vel_folder + y_list[ind])
        #v = rioxr.open_rasterio(vel_folder + y_list[ind][:-5] + 'mag' + y_list[ind][-4:])
        
        # NOTE: this doesnt work properly??
        if crop_vel:
            # Read shapefile and ensure CRS matches
            cropbox = gpd.read_file('/Volumes/Sandisk4TB/PhD_MS/TARSAN/QGIS/dagger_polygon.shp')
            cropbox = cropbox.to_crs(vx.rio.crs)
            
            # Crop using shapefile geometry
            vx = vx.rio.clip(cropbox.geometry, cropbox.crs, drop=True)
            vy = vy.rio.clip(cropbox.geometry, cropbox.crs, drop=True)
            #v = v.rio.clip(cropbox.geometry, cropbox.crs, drop=True)

    
        # get lat and lon coordinate
        y = vx.y
        x= vx.x
        
        # convert to numpy 
        vx = vx.to_numpy()
        vx = vx[0]*365     # convert from m/day to m/yr
        vy = vy.to_numpy()
        vy = vy[0]*365
       # v = v.to_numpy()
        #v = v[0]*365
        


        # load the image:
        
        # open the tiff as a raster to get bounds for plotting extent
        
        # use the base image as extent 
        rast = rasterio.open(img_folder + img_file)
        
        # extent for plotting
        extent = [rast.bounds[0], rast.bounds[2], rast.bounds[1], rast.bounds[3]]
        
        if gash == 1:
            gash_s1 = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Gash_Sentinel1/s1a-iw-grd-hh-20230108t043640-20230108t043705-046687-05989a-001_ps_ROI_gash.tiff'
            gash_rast = rasterio.open(gash_s1)
            gash_extent = [gash_rast.bounds[0], gash_rast.bounds[2], gash_rast.bounds[1], gash_rast.bounds[3]]
        
        if load_img:        
            img = plt.imread(img_folder + img_file)
            print('loaded image')

        # Calculate the 2D mean of each time slice falling in the same year
        
        if principal_directions:
            red = 10
            #red = 7
            #red = 3
            
        elif divergence:
            red = 3
        
        # slice from start to end in "red" number of chunks
        # vx_red = vx[::red,::red].to_numpy()
        # vy_red = vy[::red,::red].to_numpy()
        # x_red = x[::red].to_numpy()
        # y_red = y[::red].to_numpy()
        
        vx_red = vx[::red,::red]
        vy_red = vy[::red,::red]
        x_red = x[::red]
        y_red = y[::red]
        
        dx = np.max(np.diff(x_red))
        dy = np.max(np.diff(y_red))

        # calculate derivatives and strain rates
        
        dvxdx = (vx_red[1:-1,2:]-vx_red[1:-1,:-2])/(2*dx)
        dvxdy = (vx_red[2:,1:-1]-vx_red[:-2,1:-1])/(2*dy)
        dvydx = (vy_red[1:-1,2:]-vy_red[1:-1,:-2])/(2*dx)
        dvydy = (vy_red[2:,1:-1]-vy_red[:-2,1:-1])/(2*dy)
        
        exx = dvxdx
        eyy = dvydy
        exy = 0.5*(dvxdy+dvydx)
        
        # incompressibility 
        ezz = -(exx + eyy)
        
        e1 = 0.5*(dvxdx+dvydy) + 0.5*np.sqrt((dvxdx-dvydy)**2+4*exy**2)
        e2 = 0.5*(dvxdx+dvydy) - 0.5*np.sqrt((dvxdx-dvydy)**2+4*exy**2)

        # principal directions
        
        if principal_directions:
            print('plotting principal directions')
        
            n1 = np.sqrt(exy**2+(dvxdx-e1)**2)
            pd1x = exy/n1
            pd1y = (dvxdx-e1)/n1
            
            n2 = np.sqrt(exy**2+(dvxdx-e2)**2)
            pd2x = exy/n2
            pd2y = (dvxdx-e2)/n2
            
            print('calculated strain rates and principal directions')
            
                       
            col1 = 'r' # colour for positive (compression)
            col2 = 'b' # colour for negative (extension)
            
            #plt.figure(figsize=(16,12))
            fig, ax = plot_strainrates(x_red, y_red, e1,e2,pd1x,pd1y,pd2x,pd2y, col1, col2)
            
        elif divergence:
            print('plotting ezz')
            
            fig, ax = plt.subplots(figsize=(8, 6))
            
            ax.imshow(img, vmin= 0, vmax = 400, cmap = 'gray', extent = extent) 
            ax.set_xlim(-1.60494e6, -1.56003e6)
            ax.set_ylim(-4.6592e5, -4.2995e5)

            # Define a symmetric range for the color scale
            
            # Create plot
            c = ax.pcolormesh(x_red[1:-1], y_red[1:-1], ezz, shading='auto', cmap='RdBu', vmin=-0.5e-2, vmax=0.5e-2, alpha=0.5)
            
            # Labels and title
            #ax.set_title('Vertical Strain Rate $\\dot{\\varepsilon}_{zz}$', fontsize=14)
            #ax.set_xlabel('x [m]')
            #ax.set_ylabel('y [m]')
            #ax.set_aspect('equal')
            
            # Add colorbar
            cb = fig.colorbar(c, ax=ax)
            cb.set_label('$\\dot{\\varepsilon}_{zz}$ [1/m]', fontsize=10)
            
            plt.tight_layout()
            plt.show()

        
        if gash == 1:
            # Zoom in to the gash_extent region
            #ax1 = plt.gca()
            ax.set_xlim(-1582000, -1571000)
            ax.set_ylim(-447000, -438000)
            ax.set_aspect('equal')
            
        elif channel == 1:
            # Zoom in 
            #ax1 = plt.gca()
            ax.set_xlim(-1590195.76485995,-1565214.19730687)
            ax.set_ylim(-438590.788246874, -427884.402152700)
            ax.set_aspect('equal')
            ax.set_title(str(date))
            plt.show()
            
        elif daggers == 1:
            # Zoom in 
            ax.set_xlim(-1.58687e6, -1.57220e6)
            ax.set_ylim(-4.5101e5, -4.3616e5)
            ax.tick_params(labelsize=10)
            ax.set_title(vel[7:32],fontsize=10)
        
        else:
            #ax1 = plt.gca()
            ax.tick_params(labelsize=10)
            ax.set_title(vel[7:32],fontsize=10)
            
            
            
        
        if savepng == 1:
        
            #fig = plt.gcf()
        
            # Save the figure as a PNG file
            if principal_directions:
                filename = f"{y_list[ind][:-5]}_principal_directions.png"
                
                if daggers:
                    filename = f"{y_list[ind][:-5]}_principal_directions_daggers.png"
                    fig.savefig(vel_folder + 'e_principal_directions_daggers/' + filename, dpi=300)
                    
                else:
                    fig.savefig(vel_folder + 'e_principal_directions/' + filename, dpi=300)
                
            elif divergence:
                filename = f"{y_list[ind][:-5]}_ezz.png"
                fig.savefig(vel_folder + filename, dpi=300)
        
            # Close the figure to free memory
            plt.close(fig)
        
            print(f"Saved: {filename}")
# %%
