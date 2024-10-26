#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:26:44 2024

@author: msharp
"""

# Script to extract S1 images from the zip folders downlaoded from the ASF
# and then crop them to an area of interest

# Output file naming: e.g. s1a-iw-grd-hh-20141010t043547-20141010t043611-002762-0031af-001_ps_ROI.tiff

import os
import zipfile
import sys

module_path = '/Users/msharp/GitHub/mytoolbox/general/'
if module_path not in sys.path:
    sys.path.append(module_path)
    
from rastertools import make_outfn, assign_crs, reproject_raster, resample_raster, crop_to_extent

# %% Set up paths

raw_dir = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Sentinel-1_202206-202410/'    # folder that holds downloaded zip files from ASF
cropped_dir = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'  # folder to hold cropped rasters

reference_rast = 's1a-iw-grd-hh-20141010t043547-20141010t043611-002762-0031af-001_ps_ROI.tiff'  # filename of cropped raster that defines the extent of the output
reference_raster_path = os.path.join(cropped_dir, reference_rast)

# %% Unzip & extract S1 image

for fn in os.listdir(raw_dir):
    
    # get the path of all the zip folders
    if fn.endswith('.zip') and not fn.startswith('._'):
        zip_path = os.path.join(raw_dir, fn)
        
        # open the zip file as read-only
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            
            # Get the path of the tiff file
            for file in zip_ref.namelist():
                if 'measurement/' in file and file.endswith('.tiff'):
                    
                    print(f'Working on {file}')
                    
                    step = ''
                    
                    filename = os.path.split(file)[1]  # the filename of the current tiff
                    
                    # Extract the tiff file to the folder path
                    zip_ref.extract(file, raw_dir)
                    
                    # Reproject, crop, resample the tiff
                    input_raster_path = os.path.join(raw_dir, file)  # path to the .tiff file
                           
                    current_rast_path = input_raster_path
                    
                    # make sure the .tiff has a crs... the ASF S1 downloads do not!
                    step += '_crs'
                    outfile = make_outfn(filename, step)
                    output_raster_path = os.path.join(cropped_dir, outfile)
                    assign_crs(input_raster_path, output_raster_path)
                    
                    current_rast_path = output_raster_path
                    
                    step += '_epsg3031'
                    outfile = make_outfn(filename, step)
                    output_raster_path = os.path.join(cropped_dir, outfile)
                    reproject_raster(current_rast_path, output_raster_path, epsg = 3031)
                    
                    os.remove(current_rast_path)
                    current_rast_path = output_raster_path
                    
                    step += '_resample'
                    outfile = make_outfn(filename, step)
                    output_raster_path = os.path.join(cropped_dir, outfile)
                    resample_raster(current_rast_path, output_raster_path, target_res = 10, method = "bilinear")
                    
                    os.remove(current_rast_path)
                    current_rast_path = output_raster_path
                    
                        
                    step = '_ROI'
                    outfile = make_outfn(filename, step)
                    output_raster_path = os.path.join(cropped_dir, outfile)
                    crop_to_extent(current_rast_path, reference_raster_path, output_raster_path)
                    
                    os.remove(current_rast_path)
                