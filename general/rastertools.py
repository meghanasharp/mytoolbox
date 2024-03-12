#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:16:30 2023

@author: msharp
"""

# %%
import os
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject

# %%
def crop_raster(input_raster, reference_raster, output_path):
    # Open the reference raster to get the target extent
    with rasterio.open(reference_raster) as ref_src:
        ref_extent = ref_src.bounds
        
        print('working on raster:' + input_raster)

    # Open the input raster to crop
    with rasterio.open(input_raster) as src:
        
        input_crs = src.crs
        ref_crs = ref_src.crs
        
        # Reproject the input raster if CRS is different
        if input_crs != ref_crs:
            
            # Reproject the input raster to match the CRS of the reference raster
                
            transform, width, height = calculate_default_transform(
                input_crs, ref_crs, src.width, src.height, *src.bounds
            )
        
            data = src.read(
                out_shape=(src.count, height, width),
                resampling=Resampling.bilinear
            )
        
            # Reproject the data
            print('Reprojecting input raster...')
            reproject(
                source=data,
                destination=rasterio.band(src,1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=ref_crs,
                resampling=Resampling.bilinear
            )
        # Calculate the window to crop
        window = from_bounds(*ref_extent, transform=src.transform)
        print('Window created from bounds')

        # Read the data within the window
        data = src.read(window=window)

        # Update the metadata with the new window
        transform = src.window_transform(window)
        profile = src.profile
        profile.update({
            'width': window.width,
            'height': window.height,
            'transform': transform
        })

        # Create a new raster file with the cropped data
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(data)
            print('Saving clipped raster file')

### Example usage
# input_directory = 'path/to/yourdirectory/'
# reference_raster_path = 'path/to/your/reference_raster.tif'
# output_directory = 'path/to/your/outputdirectory/.tif'

# # List all TIFF files in the input directory
# input_raster_files = [f for f in os.listdir(input_directory) if f.endswith('.tif')]

# # Crop each raster file
# for file in input_raster_files:
#     input_raster_path = os.path.join(input_directory, file)
    
#     outfile = file[:-4] + '_cropped' + file[-4:]
#     output_raster_path = os.path.join(output_directory, outfile)

#     crop_raster(input_raster_path, reference_raster_path, output_raster_path)
    

# %% Apply to Thwaites S1 images

# input_directory = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Sentinel-1/'
# reference_raster_path = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/s1a-iw-grd-hh-20211102t043637-20211102t043702-040387-04c983-001_ps_ROI.tiff'
# output_directory = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'

# # List all TIFF files in the input directory
# input_raster_files = [f for f in os.listdir(input_directory) if f.endswith('.tiff')]

# # Crop each raster file
# for file in input_raster_files:
#     input_raster_path = os.path.join(input_directory, file)
    
#     outfile = file[:-5] + '_cropped' + file[-5:]
#     output_raster_path = os.path.join(output_directory, outfile)

#     crop_raster(input_raster_path, reference_raster_path, output_raster_path)