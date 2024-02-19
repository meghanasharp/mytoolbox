#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:52:01 2024

@author: msharp
"""

import rasterio
from rasterio.enums import Resampling

def resample_img(input_tif, output_tif, upscale_factor):
# Downsampling to 1/2 of the resolution can be done with upscale_factor = 1/2
# Resample method current hard coded in but can be changed with
# options chosen from here: https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling

    with rasterio.open(input_tif, mode="r+") as dataset:
    
        # resample data to target shape
        data = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling= Resampling.bilinear
        )
    
        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / data.shape[-1]),
            (dataset.height / data.shape[-2])
        )
        
        # Update the transform in the dataset
        dataset.transform = transform

        # Write the resampled data with updated transform to a new raster file
        with rasterio.open(output_tif, 'w', driver='GTiff', height=data.shape[1], width=data.shape[2], count=dataset.count, dtype=data.dtype, crs=dataset.crs, transform=transform) as dst:
    
           dst.write(data)
           
    return output_tif
