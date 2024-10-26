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

def assign_crs(input_raster, output_raster):
    
    print('No CRS found. Reprojecting to add CRS...') 
    
    # Open the raster file
    with rasterio.open(input_raster) as src:
        
        # Get GCPs
        gcps, gcp_crs = src.gcps
        
        if gcps:
            # Create a CRS object (replace with the correct EPSG code)
            new_crs = gcp_crs  
    
            # Create a transform from the GCPs
            transform = from_gcps(gcps)
    
            # Create a new raster with the updated CRS and transform
            metadata = src.meta.copy()
            metadata.update({
                'crs': new_crs,
                'transform': transform
            })
    
            with rasterio.open(output_raster, 'w', **metadata) as dst:
                dst.write(src.read(1), 1)  # Assuming you want to write the first band

    
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
        with gdal.Open(ref_path) as ref:

            # Get the bounding box of the reference raster
            ref_gt = ref.GetGeoTransform()
            xmin = ref_gt[0]
            ymax = ref_gt[3]
            xmax = xmin + ref_gt[1] * ref.RasterXSize
            ymin = ymax + ref_gt[5] * ref.RasterYSize
            
            bounds = (xmin, ymin, xmax, ymax)
            
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

# input_directory = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Sentinel-1/'
# reference_raster_path = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/s1a-iw-grd-hh-20211102t043637-20211102t043702-040387-04c983-001_ps_ROI.tiff'
# output_directory = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Imagery/Clipped_Data/'

# # List all TIFF files in the input directory
# input_raster_files = [f for f in os.listdir(input_directory) if f.endswith('.tiff')]

# def make_outfn(fn, step):
#     return fn[:-5] + step + fn[-5:]

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
    
#     os.remove(current_rast_path)
    

# %% GRAVEYARD

# def reproject_raster(input_raster, ref_crs, outpath):
    
#     with rasterio.open(input_raster) as src:
    
#         print(f'Reprojecting from {src.crs} to {ref_crs}')
        
        
#         # Reproject the input raster to match the CRS of the reference raster
#         transform, width, height = calculate_default_transform(
#             src.crs, ref_crs, src.width, src.height, *src.bounds
#         )
        
#         # Define the destination profile with the new CRS
#         kwargs = src.meta.copy()
#         kwargs.update({
#             'crs': ref_crs,
#             'transform': transform,
#             'width': int(width),
#             'height': int(height)
#         })
    
#         with rasterio.open(outpath, 'w', **kwargs) as dst:
                                      
#             reproject(
#                 source=rasterio.band(src,1),
#                 destination=rasterio.band(dst,1),
#                 src_transform=src.transform,
#                 src_crs=src.crs,
#                 dst_transform=transform,
#                 dst_crs=ref_crs,
#                 resampling=Resampling.bilinear
#             )


# def resample_raster(raster_path, target_resolution, resampling_method=Resampling.bilinear):
#     """
#     Resamples a raster to a specified resolution.

#     Parameters:
#     - input_raster_path (str): Path to the input raster file.
#     - output_raster_path (str): Path to save the resampled raster.
#     - target_resolution (float): Target resolution for the raster (in units of the CRS).
#     - resampling_method (Resampling): Method of resampling (default is bilinear).
#     """
    
#     # Open the original raster to read
#     with rasterio.open(raster_path) as src:
#         # Calculate the new shape based on the target resolution
#         scale_x = src.res[0] / target_resolution
#         scale_y = src.res[1] / target_resolution
#         new_width = int(src.width * scale_x)
#         new_height = int(src.height * scale_y)

#         # Update metadata for the new shape and transform
#         metadata = src.meta.copy()
#         metadata.update({
#             'width': new_width,
#             'height': new_height,
#             'transform': src.transform * src.transform.scale(scale_x, scale_y)
#         })

#         # Resample the data for each band
#         resampled_data = []
#         for i in range(1, src.count + 1):
#             data = src.read(i, out_shape=(new_height, new_width), resampling=resampling_method)
#             resampled_data.append(data)

#     # Write the resampled data back to the original file
#     with rasterio.open(raster_path, 'w', **metadata) as dst:
#         for i, data in enumerate(resampled_data, start=1):
#             dst.write(data, i)

#     print(f"Raster at {raster_path} has been resampled and overwritten with the new resolution.")

# def crop_resample_raster(input_raster, reference_raster, output_path, target_resolution = 10):
#     # set target_resolution = [] if you don't want to resample. In units of the CRS 
    
#     print('Working on raster:' + input_raster)
    
#     current_raster = input_raster
    
#     # Open the reference raster to get the target extent and reference coordinate system
#     with rasterio.open(reference_raster) as ref_src:
#         ref_extent = ref_src.bounds
#         ref_crs = ref_src.crs
    
#     # Make sure the raster has a crs 
#     with rasterio.open(current_raster, 'r+') as src:
#         if not src.crs:
#             assign_crs(current_raster, output_path)
            
#             # Now, work on the output path
#             current_raster = output_path
                           
#     # Open the raster that has the crs assigned and reproject it to the reference raster
#     with rasterio.open(current_raster, 'r+') as src:
            
#         input_crs = src.crs  # get the crs
        
#         # Reproject the input raster if CRS is different
#         if input_crs != ref_crs:      
#             reproject_raster(current_raster, ref_crs, output_path)
            
#             current_raster = output_path
                       
#     # Then, open the reprojected raster and resample it
#     if target_resolution:
#         with rasterio.open(current_raster, 'r+') as src:
#             resample_raster(current_raster, target_resolution)
            
            
#     # Finally, open the raster that is now in the correct coordinate system and crop it!
#     with rasterio.open(output_path, 'r+') as src:
        
#         # Make sure the reprojection worked before proceeding
#         print(f'Input CRS: {src.crs}, Ref CRS: {ref_crs}')
               
#         # Calculate the window to crop
#         window = from_bounds(*ref_extent, transform=src.transform)
#         print('Window created from bounds')
        
#         # Read the data within the window
#         data = src.read(window=window)
        
#         # Update the metadata with the new window
#         transform = src.window_transform(window)
#         profile = src.profile
#         profile.update({
#             'width': window.width,
#             'height': window.height,
#             'transform': transform
#         })
        
#         # Create a new raster file with the cropped data
#         with rasterio.open(output_path, 'w', **profile) as dst2:
#             print('Saving clipped raster file')
            
#             # Write the data in the window to the new file
#             try:
#                 dst2.write(data)
#                 print(f'Successfully saved cropped raster to {output_path}')
#             except:
#                 print('Error saving cropped raster... possibly no overlap?')


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
    

