#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:40:27 2023

@author: msharp
"""

import os
import shutil
import re

#WEEK 4: CHALLENGE QUESTION
#Option 3. 

# Define a function that accepts these arguments:
#     path name for source directory ("source", type: str)
#     path name for destination directory ("dest", type: str)
#     a search term string ("string", type: str)

#Example:
# source = "/Volumes/Extreme_SSD/PhD_MS/Courses/GEOG 562 Programming for Geospatial Analysis/"
# dest = "/Volumes/Extreme_SSD/PhD_MS/Courses/GEOG 562 Programming for Geospatial Analysis/test/"
# string = "test"   '
# search_and_copy(source, dest, string)


# The function finds all the files in the source directory
#    that have the “.txt” extension, reads them, and 
#    determines which ones have the desired search
#    string in their text.  
#   It then copies those files to the destination directory
#     and finally returns a list of all of the file names 
#    it copied. 

# Define a function with the inputs as above
def search_and_copy(source, dest, string):
    
    #cheeck if destination directory already exists. If not, make one!
    if os.path.exists(dest) == False:
        os.mkdir(dest)
        
    #create a list of all files and directories in the source directory
    dir_list = os.listdir(source)
    
    #initate a list to save the copied filenames in
    copiedfiles = []
    
    #loop through each file in the list of files and directories
    for file in dir_list:
        
        #if the file ends with ".txt"...
        if file.endswith(".txt"):
        
            #open file in read mode
            txtfile = open(source + file, "r") 
            
            #read the content of the file
            content = txtfile.read()
            
            #check if the string is present in the file
            if string in content:
                
                #copy the file to destination directory and keep the filename the same
                #syntax: shutil.copy(source, destination)
                shutil.copy(source + file, dest + file)
                
                #save the filename in a list
                copiedfiles += file +", "
            
            
    #return the list of filenames that were copied
    return copiedfiles
        