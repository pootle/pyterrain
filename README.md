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

Unfortunately python support for nested zips is a bit hairy, so you need to unzip the main zip file into the contained folder structure, each tile is then itself a zipfile with a few files in it. The main data is still zipped so it is small and fast to access. Use this command (or similar)

> unzip terr50_gagg_gb.zip

The code below assumes the data files are in the folder 'ost50grid'. They can be in subfolders as from the unzip above, the setup will find them.

Installation
============

There is a single python module that does all the work. 

Usage
=====

With python3:
> import pyterrain
> osd=pyterrain.osTerraintileindex('ost50grid')
> osd.summary()

This returns a brief report on what has been loaded

> osd.makeSTLfile(corner1=(743623,204635), corner2=(762823,226315), northfirst:True, fileout='glencoe.stl', 'step':1, 'mode':'bin')

The creates a binary stl file of the Glencoe area, from the westerm edge of Rannoch moor to part of Loch Leven. This file can be loaded into slic3r, or viewed with blender or other stl rendering programs.

This one does the whole UK, at max reolution. The file is big and it takes a while (12-15 mins) and generates a 13Gb file.

> osd.smartSTLfile(corner1=(0,0), corner2=(1230000,660000), northfirst=True, fileout='all.stl'

You can derive corner values from streetmap.co.uk, or use lat long from other maps, like openstreetmap, google maps, etc.

Use help(osd.smartSTLfile) for details of the parameters. For example, if you want a slightly more human readable file, use mode='txt as a parameter, this is slightly slower and creates a much larger text file.
