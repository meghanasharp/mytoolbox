#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:16:30 2023

@author: msharp
"""

# %%
import os
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
from rasterio.transform import from_gcps

from osgeo import gdal, osr

# %%

# def assign_crs(input_raster, output_raster):
    
#     print('No CRS found. Reprojecting to add CRS...') 
    
#     # Open the raster file
#     with rasterio.open(input_raster) as src:
        
#         # Get GCPs
#         gcps, gcp_crs = src.gcps
        
#         if gcps:
#             # Create a CRS object (replace with the correct EPSG code)
#             new_crs = gcp_crs  
    
#             # Create a transform from the GCPs
#             transform = from_gcps(gcps)
    
#             # Create a new raster with the updated CRS and transform
#             metadata = src.meta.copy()
#             metadata.update({
#                 'crs': new_crs,
#                 'transform': transform
#             })
    
#             with rasterio.open(output_raster, 'w', **metadata) as dst:
#                 dst.write(src.read(1), 1)  # Assuming you want to write the first band


def assign_crs(input_raster, output_raster):
    
    print('No CRS found. Reprojecting to add CRS...') 
    
    # Open the raster file
    with rasterio.open(input_raster) as src:
        
        # Get GCPs
        gcps, gcp_crs = src.gcps
    
    with gdal.Open(input_raster) as src:
        if gcps:
            
            output_rast = gdal.Warp(output_raster, src, dstSRS=gcp_crs)
            
            # Close the file
            output_rast = None

    
def reproject_raster(input_path, output_path, epsg = 3031):
    
    print(f'Reprojecting raster to EPSG:{epsg}')
    
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(epsg)
    
    with gdal.Open(input_path) as src:       
        output_rast = gdal.Warp(output_path, src, dstSRS=target_srs)
        
        # Close the file
        output_rast = None
        
        
def resample_raster(input_path, output_path, target_res = 10, method = "bilinear"):
    
    print(f'Resampling raster to {target_res}')

    with gdal.Open(input_path) as src:     
        
        # Resample the raster
        output_rast = gdal.Warp(
            output_path,
            src,
            xRes=target_res,  # Target resolution in X direction
            yRes=target_res,  # Target resolution in Y direction
            resampleAlg=method  # Resampling method
        )
        
        # Close the file
        output_rast = None
        
def crop_to_extent(input_path, ref_path, output_path):
    
    print('Cropping raster...')
    
    # Open both rasters
    with gdal.Open(input_path) as src:
        with rasterio.open(ref_path) as ref:
            
            bounds = ref.bounds
            
            # Crop the input raster to match the target's extent
            output_raster = gdal.Warp(
                output_path,
                src,
                outputBounds= bounds
            )
            
            # Close the file
            output_raster = None
            
            output_fn = os.path.split(output_path)[1]
            print(f'Cropped raster to: {output_fn}')

def make_outfn(fn, step):
    return fn[:-5] + step + fn[-5:]

# %% Apply to Thwaites S1 images

# input_directory = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/test/'
# reference_raster_path = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/s1a-iw-grd-hh-20211102t043637-20211102t043702-040387-04c983-001_ps_ROI.tiff'
# output_directory = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/test_clipped/'

# # List all TIFF files in the input directory
# input_raster_files = [f for f in os.listdir(input_directory) if f.endswith('.tiff') and not f.startswith('._')]

# step = ''
# # Crop each raster file
# for file in input_raster_files:
#     input_raster_path = os.path.join(input_directory, file)
    
#     current_rast_path = input_raster_path
    
#     step += '_crs'
#     outfile = make_outfn(file, step)
#     output_raster_path = os.path.join(output_directory, outfile)
#     assign_crs(input_raster_path, output_raster_path)
    
#     current_rast_path = output_raster_path
    
#     step += '_epsg3031'
#     outfile = make_outfn(file, step)
#     output_raster_path = os.path.join(output_directory, outfile)
#     reproject_raster(current_rast_path, output_raster_path, epsg = 3031)
    
#     os.remove(current_rast_path)
#     current_rast_path = output_raster_path
    
#     step += '_resample'
#     outfile = make_outfn(file, step)
#     output_raster_path = os.path.join(output_directory, outfile)
#     resample_raster(current_rast_path, output_raster_path, target_res = 10, method = "bilinear")
    
#     os.remove(current_rast_path)
#     current_rast_path = output_raster_path
    
#     step = '_ROI'
#     outfile = make_outfn(file, step)
#     output_raster_path = os.path.join(output_directory, outfile)
#     crop_to_extent(current_rast_path, reference_raster_path, output_raster_path)
    

#     #os.remove(current_rast_path)


