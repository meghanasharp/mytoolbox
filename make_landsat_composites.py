#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:17:32 2024

@author: msharp
"""

from osgeo import gdal
import os
import resample

### Meghan Sharp Jan 16/ 2024
### This script will create composite images for all Landsat7 and Landsat8 data in a specified folder

# Path of directory containing Landsat Images
path = '/Volumes/Sandisk4TB/PhD_MS/SCAR/Landsat_Imagery/new/'

#root = 'LE07_L2SR_218106_20130222_20200908_02_T2_SR'

# Get a list of files in the directory
file_list = os.listdir(path)

# Filter files that end with 'MTL.txt'
mtl_files = [file for file in file_list if file.endswith('MTL.txt') and not file.startswith('._')]

for img in mtl_files:
     
    root = img[:-8]

    satellite = root[0:4]

    if satellite == 'LE07':
        redband= path + root + '_SR_B3.TIF'
        greenband = path + root + '_SR_B2.TIF'
        blueband = path + root + '_SR_B1.TIF'
        
        print('Preparing composite: ', redband, ', ', greenband, ', ', blueband)
        
        output_vrt = path + 'Landsat_Composites/' + root + '_SR_B3_B2_B1.VRT'
    
    elif satellite == 'LC08':
        redband= path + root + '_SR_B4.TIF'
        greenband = path + root + '_SR_B3.TIF'
        blueband = path + root + '_SR_B2.TIF'
        
        output_vrt = path + 'Landsat_Composites/' + root + '_SR_B4_B3_B2.VRT'
    
    output_tif = output_vrt[:-4] + '.TIF'

    
    # Create a virtual raster (VRT) with option "separate=True" to stack the images as separate bands
    print('Creating virtual raster...')
    my_vrt = gdal.BuildVRT(output_vrt, [redband, greenband, blueband], separate=True, callback = gdal.TermProgress_nocb)
    
    # Translate virtual raster (VRT) into GeoTiff
    print('Translating VRT to GeoTiff...')

    gdal.Translate(output_tif, my_vrt, format='GTiff',
                       callback=gdal.TermProgress_nocb)
    
    #Downsample the compsoite raster for faster processing in QGIS, etc
    print('Downsampling...')
    resampled_tif = output_tif[:-4] + '_downsampled' + '.TIF'
    resample.resample_img(output_tif, resampled_tif, 1/2)
    
    print('Downsampled, composite raster save to: ', resampled_tif)

