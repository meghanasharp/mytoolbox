#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:41:03 2024

@author: msharp
"""
import matplotlib.pyplot as plt
import rasterio


def nice_quiver(lon, lat, vx, vy, v, vmin, vmax, units = 'm/yr', n =30, img = '', title = '', savefig = 0, outname = 'fig.png', cmap = 'magma_r'):

    ##############
    # This function was built with the intention of itslive data, 
    # but will probably work for any velocity grid.
    # 
    # INPUTS:
        # lon: longitude coordinates as a meshgrid (numpy array)
        # lat: longitude coordinates as a meshgrid (numpy array)
        # vx: x-component of velocity with same grid dimesions as lat and lon (numpy array)
        # vy: y-component of velocity with same grid dimesions as lat and lon (numpy array)
        # v: magntidue of velocity with same grid dimesions as lat and lon (numpy array)
        # vmin: minimum value for velocity magnitude colours (float)
        # vmax: maxmimum value for velocity magnitude colours (float)
        # n: spacing of arrows (int, DEFAULT is every 30 grid cells)
        # img: path to geotiff for base image. (str, OPTIONAL)
        # title: title for figure (str, OPTIONAL)
        # savefig: set = 1 if you want to save the figure (bool, OPTIONAL)
        # outname: output filename for figure to save to (str, DEFAULT is fig.png)
        # cmap: colourmap code for velocity magnitude (str, DEFAULT is 'magma_r')
    #
    # OUTPUTS:
        # returns the figure handle
    #
    # Author: Meghan A. Sharp, Feb/ 2024
    #############
    
    # initiate the figure
    f, ax = plt.subplots()
    
    # open the base image
    img_folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'
    img_file = 's1a-iw-grd-hh-20210530t043630-20210530t043655-038112-047f7a-001_ps_ROI.tiff'
    img = plt.imread(img_folder + img_file)
    
    # open the tiff as a raster to get bounds for plotting extent
    rast = rasterio.open(img_folder + img_file)
    extent = [rast.bounds[0], rast.bounds[2], rast.bounds[1], rast.bounds[3]]
    
    plt.imshow(img, vmin= 0, vmax = 400, cmap = 'gray', extent = extent)
    vplot = plt.imshow(v, vmin= vmin, vmax = vmax, cmap = cmap, extent=[lon.min(), lon.max(), lat.min(),lat.max()], alpha = 0.6)
    plt.quiver(lon[::n, ::n], lat[::n, ::n], vx[::n, ::n], vy[::n, ::n], scale = 0.3, units = 'xy')  # python slicing is array[start:stop:step]
    
    
    # set figure properties 
    plt.title(title)
    plt.xlim([rast.bounds[0], rast.bounds[2]])
    plt.ylim([rast.bounds[1], rast.bounds[3]])
    ax.set_xticks([])
    ax.set_yticks([])
    
    cbar = plt.colorbar(vplot)
    cbar.set_label('v magnitude ' + units)
    #plt.show()
    
    # save the figure
    if savefig:
        plt.savefig(outname, dpi = 150)
        
    return f