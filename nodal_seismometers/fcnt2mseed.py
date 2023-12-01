# -*- coding: utf-8 -*-

"""
Standalone script to convert raw Fairfield fcnt files to mseed files.

For more information on how to read fcnt files with the ObsPy library, check
following link: https://docs.obspy.org/master/packages/obspy.io.rg16.html

For more information on how to write mseed files with the ObsPy library, check
following link: https://docs.obspy.org/packages/obspy.io.mseed.html

Maeva Pourpoint - IRIS/PASSCAL
"""

# Imports
from obspy.core import read
from pathlib import Path

# Change path to fcnt files if necessary
path2fcnt = Path('/Volumes/Extreme_SSD/PhD_MS/TARSAN/Nodes/RAW-NODE')
#fcnt_files = path2fcnt.glob('*.fcnt')
fcnt_files = path2fcnt.glob('1.0.0.fcnt')

# Create output directory (Change directory if necessary)
pathout = Path('/Volumes/Extreme_SSD/PhD_MS/TARSAN/Nodes/NODE-MINISEED')
if not pathout.is_dir():
    pathout.mkdir()

# Loop over fcnt files
for fcnt_file in fcnt_files:
    print(f"Converting fcnt file - {fcnt_file.name} - to a mseed format")
    # Read each fcnt file
    st = read(str(fcnt_file))
    # Convert raw data to mseed format
    fileout = pathout.joinpath(fcnt_file.name.replace('fcnt', 'mseed'))
    st.write(fileout, format='MSEED')
