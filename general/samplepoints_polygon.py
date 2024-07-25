#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 09:24:59 2024

@author: msharp
"""

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import transform
from pyproj import Proj, Transformer
import pandas as pd

# %% Function to sample points along the polygon boundary at specified interval
def sample_points(polygon, interval=0.1):
    points = []
    distance = 0.0

    while distance < polygon.length:
        point = polygon.interpolate(distance)
        points.append((point.x, point.y))
        distance += interval

    return points

# %% Read the shapefile
shapefile_path = '/Volumes/Sandisk4TB/PhD_MS/LeConte/sprieberg_outlines/spireberg_outlines.shp'
gdf = gpd.read_file(shapefile_path)

# Initialize a list to hold the sampled points with their coordinates
data = []

# Define the UTM zone 8N projection
utm_proj = Proj(proj='utm', zone=8, ellps='WGS84', south=False)

# Define the UTM zone 8N transformer
transformer_utm = Transformer.from_crs("epsg:3857", "epsg:32608", always_xy=True)
transformer_wgs84 = Transformer.from_crs("epsg:3857", "epsg:4326", always_xy=True)

# %% Iterate through each polygon in the shapefile
for index, row in gdf.iterrows():
    polygon_id = row['id']  
    polygon = row.geometry
    
    # Sample points along the boundary
    sampled_points = sample_points(polygon.boundary)
    
    # Initialize a list to hold the sampled points with their coordinates for this polygon
    data = []
    
    for x, y in sampled_points:
        # Convert to UTM
        utm_x, utm_y = transformer_utm.transform(x, y)
        
        lon, lat = transformer_wgs84.transform(x, y)
        
        # Append the coordinates to the data list
        data.append([lat, lon, utm_x, utm_y])
    
    # Create a DataFrame from the data
    df = pd.DataFrame(data, columns=['Latitude', 'Longitude', 'UTM_X', 'UTM_Y'])
    
    # Write the DataFrame to a CSV file
    csv_output_path = '/Volumes/Sandisk4TB/PhD_MS/LeConte/sprieberg_outlines/' + f'{polygon_id}.csv'
    df.to_csv(csv_output_path, index=False)

    print(f"CSV file for polygon {polygon_id} has been written to {csv_output_path}")
