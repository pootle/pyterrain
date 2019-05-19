# pyterrain
Python to convert digital elevation info to stl for 3d printing.

This python reads the digital elevation data files from UK Ordnance Survey and can produce print ready stl for anywhere in the uk.

I originally used this web site: http://jthatch.com/Terrain2STL/
but I wasn't too impressed with the results, so I decided to roll my own.

I have also tried using the SRTM data as scanned from the ISS and Aster data as source data, but the most recent versions of both
these data sets (SRT V3 - 1 arc second and aster global DZEM V2) are noisy and have many anomolies (like big holes in lakes and
turning cliffs into slopes). The OS data (50m grid spacing - available for free with some license restrictions) in comparison is
far cleaner. It also has the advantage of being ready projected onto a cartesian plane so is much faster to process.

Pre-Installation
================

The code uses pyproj to map from lat/long coordinates to ordnance survey coordinates:

on linux: sudo apt install python3-pyproj.

You also need to grab the Ordnance Survey terrain 50 data: https://www.ordnancesurvey.co.uk/business-and-government/products/terrain-50.html

The code below assumes the data files are in the folder 'ost50grid'. They can be in subfolders, the setup will find them.

Don't unzip the OS zip files, the code reads directly out of the zip file.

Installation
============

There is a single python module that does all the work. 

Usage
=====

With python3:
> import pyterrain
> osd=pyterrain.osTerraintileindex('ost50grid')
> osd.summary()   # returns a brief report on what has been loaded
> osd.makeSTLfile(corner1=(743623,204635), corner2=(762823,226315), northfirst:True, fileout='glencoe.stl', 'step':1, 'mode':'bin')

This one does the whole UK, at max reolution. The file is big and it takes a while (12-15 mins) and generates a 13Gb file.

You can derive corner values from streetmap.co.uk, or use lat long from other maps, like openstreetmap, google maps, etc.
