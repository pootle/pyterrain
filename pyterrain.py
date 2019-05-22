import pathlib
import zipfile
import xml.etree.ElementTree as xxml
from collections import OrderedDict
import numpy
import pyproj
import time
import struct

def _xmlgetchild(node, cname):
    """
    a simplified fetcher that ignores namespaces

    node is an xml node (as in getroot) or any other node

    cname is a string with the name of the desired node
            - or -
             an array of strings to use in turn to find the desired node
    """
    if isinstance(cname,str):
        last=True
        tagstr=cname
    else:
        last = len(cname)==1
        tagstr=cname[0]
    tagsplit=tagstr.split(':')
    if len(tagsplit) > 1:
        tagstr=tagsplit[0]
        tagidx=int(tagsplit[1])
    else:
        tagidx=0
    curidx=0
    for n in node:
        tagend=n.tag.split('}')
        tagr=tagend[0] if len(tagend)==1 else tagend[1]
        if tagr==tagstr:
            if curidx==tagidx:
                if last:
                    return n
                else:
                    return _xmlgetchild(n,cname[1:])
            else:
                curidx+=1

class binarywriter():
    """
    write a binary stl file with normals always set to 0,0,0, and no extras
    """
    def __init__(self, file, header, trianglecount=2):
        assert trianglecount %2==0
        self.structdef='<12fH'
        self.opf = file.open(mode='wb')
        self.opf.write(struct.pack('<80s', '{:80s}'.format(header).encode()))
        self.sizepos=self.opf.tell()
        self.opf.write(struct.pack('<I',0))
        self.unitsize=struct.calcsize(self.structdef)
        self.buffsize=self.unitsize*trianglecount
        self.buff=bytearray(self.buffsize)
        for i in range(int(trianglecount/2)):
            struct.pack_into('<3f', self.buff, i*self.unitsize, 0,0,0)
            struct.pack_into('<H', self.buff, 48+i*self.unitsize, 0)
        self.buffpos=0
        self.tricount=0
        self.bufftri=trianglecount

    def writetriangle(self, *floats):
        """
        writes a triangle ( 3 vertices) from the 9,
        values given. Normal convention with outside
        faceon the anticlockwise surface
        """
        assert len(floats)==9
        self.tricount+=1
        struct.pack_into('<9f',self.buff,self.buffpos+12, *floats)
        self.buffpos+=self.unitsize
        if self.buffpos == self.buffsize:
            self.opf.write(self.buff)
            self.buffpos=0
        
    def writepaira(self, vx0, vx1, vy0, vy1, xy00, xy01, xy10, xy11):
        """
        given a quad with points on a regular grid, write out the 2 triangles.
        
        The 4 points are expanded from the given values:
        
        p1: x0, y0, xy00
        p2: x1, y0, xy01
        p3: x1, y1, xy11
        p4: x0, y1, xy10
        """
        struct.pack_into('<9f', self.buff, self.buffpos+12,
                vx0, vy0, xy00,
                vx1, vy0, xy01,
                vx1, vy1, xy11,
                )
        self.buffpos+=self.unitsize
        struct.pack_into('<9f', self.buff, self.buffpos+12,
                vx0, vy0, xy00,
                vx1, vy1, xy11,
                vx0, vy1, xy10
                )
        self.tricount+=2
        self.buffpos+=self.unitsize
        if self.buffpos >= self.buffsize:
            self.opf.write(self.buff)
            self.buffpos=0

    def close(self):
        if self.buffpos > 0:
            self.opf.write(self.buff[0:self.buffpos])
        print('file has %d triangles, size is %d bytes' % (self.tricount, self.opf.tell()))
        self.opf.seek(self.sizepos)
        self.opf.write(struct.pack('<I',self.tricount))
        self.opf.close()
        self.opf=None

class osTerraintileindex():
    """
    A container to gather info about a load of tiles.
    
    British National Grid proj=tmerc ellps=airy lat_0=49dN lon_0=2dW k_0=0.9996012717 x_0=400000 y_0=-100000
    from here https://github.com/OrdnanceSurvey/proj-4/blob/master/nad/world
    
    using this projection on 'kidsgrove':
        [-2.240827],[53.086212] => [383871.7890803295], [354393.0446889817]
        streetmap's convertion  =>  383967, 354360    (error is -96, 33)
    and on 'inverness':
        [-4.222802],[57.478378] => [266733.87969994405], [845280.3901484823]
        streetmap's convertion  => [266816, 845306]    (error is -82, -26)
    roundabout in inverness:
        google's loc 57.4919041,-4.2171973
        streetmap: 57.491758, -4.215071
        streetmap os coords:     [267328 846780
        my proj from streetmap:  [267246], [846754]   (-82,  -26)
        my proj from google:     [267119], [846775]   (-209, -5)
    """
    def __init__(self, tilefolder):
        projparams={'proj':'tmerc', 'ellps': 'airy', 
                    'lat_0': '49dN', 'lon_0': '2dW',
                    'k_0':0.9996012717, 'x_0':400000, 'y_0':-100000}
        self.proj =pyproj.Proj(**projparams)
        self.tiles={}
        d=pathlib.Path(tilefolder).expanduser()
        if not d.exists():
            raise ValueError('>%s< cannot find folder for tiles' % str(d))
        tempgrid=self.loadFolderData(d)
        northvals=list(tempgrid.keys())
        minn=min(northvals)
        maxn=max(northvals)
        anorth=tempgrid[northvals[0]]
        atile=anorth[list(anorth.keys())[0]]
        self.tileStepN=atile.ns[1]-atile.ns[0]
        self.tileStepE=atile.ew[1]-atile.ew[0]
        self.gstepN, self.gstepE=atile.getStepUnits()
        self.tileRows, self.tileCols=atile.getRowColCount()
        emptyn={}
        self.tgrid=OrderedDict(((k, tempgrid[k]) if k in tempgrid else (k, emptyn) for k in range(minn, maxn+1, self.tileStepN)))
        for rkey, rentries in self.tgrid.items():
            kl=sorted(list(rentries.keys()))
            self.tgrid[rkey]=OrderedDict(((k, rentries[k]) for k in kl))
        tlist=self.tiles.values()
        self.limS=min([t.ns[0] for t in tlist])
        self.limN=max([t.ns[1] for t in tlist])
        self.limE=min([t.ew[0] for t in tlist])
        self.limW=max([t.ew[1] for t in tlist])


    def smartSTLfile(self, corner1, corner2, fileout, useLatLon=False, northfirst=False, 
                    step=1, scale=.1, vscale=1.2, mode='bin', tribuffer=128):
        """
        write an stl file for the area from corner1 to corner2.
        
        corner1   : coordinates of 1 corner of the area 
        
        corner2   : coordinates of the opposite corner of the area
        
        fileout   : full name of the stl file to write
        
        useLatLon : if True, corners are in degrees latitude and longitude, otherwise they are OS Northing / Easting

        northfirst: order of coordinates in corners, defaults to east first
        
        step      : Specifies interval used (both lat and lon) for points processed. 
                    1 => every height point used
                    2 => every other point used (1/4 number of triangles and ~1/4 file size)
                    3 => every 3rd.....   useful when plotting large areas to avoid overly detailed meshes
        scale     : scales the output in x, y and z by the given number. The STL file has no unit of size
                    as such, although the values are in meters.
        
        vscale    : additional scaling applied to height values. 1 means heights are correct, 1.5 means they
                    are 50% higher. 

        mode      : selects stl generator:
                    'txt' : basic text stl - large file, slow to generate, readable output
                    'bin' : small file fast generation and not easily human readable 

        tribuffer : number of triangles assembled in buffer before written to file (must be even) 
        """
        clockstart=time.time()
        cpustart=time.clock()
        assert len(corner1)==2
        assert len(corner2)==2
                # first sort the start and end values so they are easting then northing
        cor1, cor2 = (list(reversed(corner1)), list(reversed(corner2))) if northfirst else (corner1, corner2)
                # convert lat / lon values to OS co-ordinates if necessary
        if useLatLon:
            xarr,yarr=self.proj([cor1[0],cor2[0]], [cor1[1], cor2[1]])
            c1=(round(xarr[0]), round(yarr[0]))
            c2=(round(xarr[1]), round(yarr[1]))
            print('input coords 1: N{:7.3f} E{:7.3f} => {:d} {:d}'.format(cor1[1],cor1[0],c1[1],c1[0]))
        else:
            c1=cor1
            c2=cor2
                # sort start and end values so end is always > start
        Estart, Eend = (c1[0], c2[0]) if c1[0] < c2[0] else (c2[0], c1[0])
        Nstart, Nend = (c1[1], c2[1]) if c1[1] < c2[1] else (c2[1], c1[1])
                # generate a list of North tile indices and get them into reverse order for the Northing list
        tilesE=list(range(int(Estart/self.tileStepE)*self.tileStepE, 
                          (int((Eend-1)/self.tileStepE)+1)*self.tileStepE, self.tileStepE))
        colstart=0 if Estart <= tilesE[0] else round((Estart-tilesE[0])/self.gstepE)
        colend  =round((Eend-tilesE[0])/self.gstepE)
        print('using columns', colstart, colend)
        tilesN=list(reversed(list(range(int(Nstart/self.tileStepN)*self.tileStepN, 
                          (int((Nend-1)/self.tileStepN)+1)*self.tileStepN, self.tileStepN))))
        print('using Ns', tilesN, 'from', Nstart, 'to', Nend)
#        lastline  = self.tileRows - (0 if Nstart<tilesN[-1] else round((Nend-tilesN[-1])/self.gstepN))
        firstline = 0 if Nend>=(tilesN[0]+self.tileStepN) else round((tilesN[0]+self.tileStepN-Nend)/self.gstepN)
        lastline  = round((tilesN[0]+self.tileStepN-Nstart)/self.gstepN)
        print('using lines', firstline, lastline)
        rowgen=self._rowmaker(tilesN, firstline, lastline, tilesE, colstart, colend, step, zscale=scale*vscale)
        xbase=self.gstepE*(colstart-colend)/2
        colxvals=[(x*self.gstepE+xbase)*scale for x in range(0, colend-colstart, step)]
        print('xvals: %d to %d step %d (%d * %d)' % (colxvals[0], colxvals[-1], step*self.gstepE, step, self.gstepE))
        if mode=='txt':
            filemaker=self.makestltextfile
        elif mode=='bin':
            filemaker=self.makestlbinfile
        else:
            raise ValueError('>%s< is not a valid mode.' % mode)
        filemaker(fileout=pathlib.Path(fileout).expanduser(),
                 rowgen=rowgen,
                 colxvals=colxvals,
                 nextyval=-(lastline-firstline)*self.gstepN/2,
                 yinc=step*self.gstepN*scale,
                 scale=scale)
        clockend=time.time()
        cpuend=time.clock()
        stlclock=clockend-clockstart
        stlcpu=cpuend-cpustart
        print('%s written in %5.2f seconds, %3.1f%% cpu' % (fileout, stlclock, stlcpu/stlclock*100))

    def makestlbinfile(self, fileout, rowgen, colxvals, nextyval, yinc, scale):
        row0=next(rowgen)
        rzl=len(row0)
        print('row0len', rzl)
        row0y=nextyval
        basey=row0y
        nextyval+=yinc
        zbase=-25
        writer=binarywriter(fileout,'Pootle terrain',trianglecount=2)
        for x in range(len(row0)-1):
            c1=(colxvals[x]  , row0y, zbase)
            c3=(colxvals[x+1], row0y, row0[x+1])
            writer.writetriangle(*c1, colxvals[x+1], row0y, zbase, *c3)
            writer.writetriangle(*c1, *c3, colxvals[x]  , row0y, row0[x])
        try:
            while True:
                row1=row0
                row1y=row0y
                row0y=nextyval
                nextyval+=yinc
                row0=next(rowgen)
                if len(row0) != rzl:
                    print('mismatched len - %d' % len(row0))
                cx1=colxvals[0]
                c1=(cx1  , row1y, zbase)
                c3=(cx1  , row0y, row0[0])
                writer.writetriangle(*c1, cx1  , row1y, row1[0], *c3)
                writer.writetriangle(*c1, *c3, cx1  , row0y, zbase)
                for x in range(len(row0)-1):
                    cx0=cx1
                    cx1=colxvals[x+1]
                    writer.writetriangle(cx0, row0y, row0[x], cx0, row1y, row1[x], cx1, row1y, row1[x+1])
                    writer.writetriangle(cx0, row0y, row0[x], cx1, row1y, row1[x+1], cx1, row0y, row0[x+1])
#                    writer.writepaira(cx0, cx1, row0y, row1y, row0[x], row1[x], row0[x+1], row1[x+1])
                c1=(colxvals[-1], row0y, zbase)
                c3=(colxvals[-1], row1y, row1[-1])
                writer.writetriangle(*c1, colxvals[-1], row0y, row0[-1], *c3)
                writer.writetriangle(*c1, *c3, colxvals[-1], row1y, zbase)
        except StopIteration:
            pass
        for x in range(len(row1)-1):
            c1=(colxvals[x]  , row1y, zbase)
            c3=(colxvals[x+1], row1y, row1[x+1])
            writer.writetriangle(*c1, colxvals[x]  , row1y, row1[x], *c3)
            writer.writetriangle(*c1, *c3, colxvals[x+1], row1y, zbase)
        c1=(colxvals[0], basey, zbase)
        c3=(colxvals[-1], row1y, zbase)
        writer.writetriangle(*c1, colxvals[0], row1y, zbase, *c3)
        writer.writetriangle(*c1, *c3, colxvals[-1], basey, zbase)
        writer.close()

    def makestltextfile(self, fileout, rowgen, colxvals, nextyval, yinc, scale):
        row0=next(rowgen)
        row0y=nextyval
        basey=row0y
        nextyval+=yinc
        zbase=-25
        tricount=0
        with fileout.open(mode='w') as pfo:
            pfo.write('solid \n')
            for x in range(len(row0)-1):
                pfo.write(facetstr1.format(
                        c1=(colxvals[x]  , row0y, zbase),
                        c2=(colxvals[x+1], row0y, zbase),
                        c3=(colxvals[x+1], row0y, row0[x+1]),
                        c4=(colxvals[x]  , row0y, row0[x]),))
                tricount+=2
            try:
                while True:
                    row1=row0
                    row1y=row0y
                    row0y=nextyval
                    nextyval+=yinc
                    row0=next(rowgen)
                    pfo.write(facetstr1.format(
                        c1=(colxvals[0]  , row1y, zbase),
                        c2=(colxvals[0]  , row1y, row1[0]),
                        c3=(colxvals[0]  , row0y, row0[0]),
                        c4=(colxvals[0]  , row0y, zbase)))
                    tricount+=2
                    for x in range(len(row0)-1):
                        pfo.write(facetstr1.format(
                                c1=(colxvals[x],   row0y, row0[x]), 
                                c2=(colxvals[x],   row1y, row1[x]),
                                c3=(colxvals[x+1], row1y, row1[x+1]), 
                                c4=(colxvals[x+1], row0y, row0[x+1])))
                        tricount+=2
                    pfo.write(facetstr1.format(
                        c1=(colxvals[-1], row0y, zbase),
                        c2=(colxvals[-1], row0y, row0[-1]),
                        c3=(colxvals[-1], row1y, row1[-1]),
                        c4=(colxvals[-1], row1y, zbase)))
                    tricount+=2
            except StopIteration:
                pass
            for x in range(len(row1)-1):
                pfo.write(facetstr1.format(
                        c1=(colxvals[x]  , row1y, zbase),
                        c2=(colxvals[x]  , row1y, row1[x]),
                        c3=(colxvals[x+1], row1y, row1[x+1]),
                        c4=(colxvals[x+1], row1y, zbase),))
                tricount+=2
            pfo.write(facetstr1.format(
                c1=(colxvals[0], basey, zbase),
                c2=(colxvals[0], row1y, zbase),
                c3=(colxvals[-1], row1y, zbase),
                c4=(colxvals[-1], basey, zbase)))
            tricount+=2
            print('%d triangles written, %d points' % (tricount, tricount*3))
            pfo.write('endsolid terrain\n')
    
    def loadFolderData(self, fpath, tempgrid={}):
        """
        Process the contents of a folder by processing each useful file.
        
        fpath: a pathlib.Path to a folder
        """
        for zfn in fpath.iterdir():
            if zfn.is_dir():
                self.loadFolderData(zfn, tempgrid=tempgrid)
            elif zfn.suffix=='.zip':
                nparts=zfn.with_suffix('').name.split('_')
                if len(nparts)==3 and nparts[1]=='OST50GRID':
                    t=ostileTerrain50(zfn)
                    self.tiles[nparts[0].upper()] = t
                    if t.ns[0] in tempgrid:
                        nlist=tempgrid[t.ns[0]]
                    else:
                        nlist={}
                    tempgrid[t.ns[0]]=nlist
                    if t.ew[0] in nlist:
                        print('duplicate location %d in file %s ignored, original file %s' % (t.ew[0], str(zfn),
                                    nlist[t.ew[0]].tilezip))
                    else:
                        nlist[t.ew[0]]=t       
        return tempgrid
        
    def _rowmaker(self, tilesN, firstline, lastline, tilesE, colstart, colend, step, zscale):
        emptyrow=None
        rowgens=None
        dummytile=b'0 '*self.tileCols
        linesleft=lastline-firstline
        for nix, nkey in enumerate(tilesN):
            skiplines=round(firstline+step*.45) if nix==0 else 0
            tilelinesleft=self.tileRows-(skiplines if nix==0 else 0) 
            if nkey in self.tgrid:
                gridrow=self.tgrid[nkey]
            else:
                gridrow={}
            rowtilereaders=[gridrow[ekey].opentile(
                    skiplines=skiplines) if ekey in gridrow else None for ekey in tilesE]
            while tilelinesleft >=0 and linesleft >=0:
                subs=[dummytile if reader is None else reader() for reader in rowtilereaders]
                subrow=b' '.join(subs).split()[colstart:colend:step]
#                    print('subrow len', len(subrow), 'tilelinesleft', tilelinesleft)
                rowints=[float(n)*zscale for n in subrow]
#                    print('first', rowints[0], 'last', rowints[-1])
                linesleft     -= 1 + step - 1
                tilelinesleft -= 1 + step - 1
                yield rowints
                for _ in range(step-1):
                    for t in rowtilereaders:
                        if not t is None:
                            t()
            for ekey in tilesE:
                if ekey in gridrow:
                    gridrow[ekey].closetile()

    def summary(self):
        return '%d tiles %d < easting < %d, %d < northing < %d, stepN: %d, stepE: %d, %d rows' % (len(self.tiles),
                            self.limE, self.limW, self.limS, self.limN, self.tileStepN, self.tileStepE, len(self.tgrid))

facetstr1=(
    'facet normal 0 0 -1\n'
    '  outer loop\n'
    '  vertex {c1[0]:6f} {c1[1]:6f} {c1[2]:6f}\n'
    '  vertex {c2[0]:6f} {c2[1]:6f} {c2[2]:6f}\n'
    '  vertex {c3[0]:6f} {c3[1]:6f} {c3[2]:6f}\n'
    '  endloop\n'
    'endfacet\n'
    'facet normal 0 0 -1\n'
    '  outer loop\n'
    '  vertex {c1[0]:6f} {c1[1]:6f} {c1[2]:6f}\n'
    '  vertex {c3[0]:6f} {c3[1]:6f} {c3[2]:6f}\n'
    '  vertex {c4[0]:6f} {c4[1]:6f} {c4[2]:6f}\n'
    '  endloop\n'
    'endfacet\n'
)

class ostileTerrain50():
    """
    a class to represent the zip file Ordnance survey uses for a Terrain50 tile
    
    Just records the bounds and size of the tile contents
    """
    def __init__(self, tilefiles):
        """
        prepares a terrain50 instance with just some basic basic info about the tile
        """
        if not zipfile.is_zipfile(str(tilefiles)):
            raise ValueError('file %s does not appear to be a zipfile')
        self.tilezip=tilefiles
        with zipfile.ZipFile(str(self.tilezip), mode='r') as zf:
            txml=xxml.parse(zf.open(tilefiles.name[0:4].upper()+'.gml')).getroot()
        env=_xmlgetchild(txml,('boundedBy', 'Envelope'))
        llc=_xmlgetchild(env,'lowerCorner').text.split(' ')
        upc=_xmlgetchild(env,'upperCorner').text.split(' ')
        self.ew=[int(llc[0]), int(upc[0])]
        self.ns=[int(llc[1]), int(upc[1])]
        self.osname=tilefiles.name[0:4]
        self.zf=None
        self.demf=None

    def getRowColCount(self):
        with zipfile.ZipFile(str(self.tilezip), mode='r') as zf:
            with zf.open(self.tilezip.name[0:4].upper()+'.gml') as xmlf:
                txml=xxml.parse(xmlf).getroot()
        ginfo=_xmlgetchild(txml, ('member', 'ElevationGridCoverage', 'rectifiedGridDomain', 'RectifiedGrid', 'limits', 'GridEnvelope'))
        counts_h=[int(n) for n in _xmlgetchild(ginfo, 'high').text.split()]
        counts_l=[int(n) for n in _xmlgetchild(ginfo, 'low').text.split()]
        return counts_h[0]-counts_l[0], counts_h[1]-counts_l[0]

    def getStepUnits(self):
        with zipfile.ZipFile(str(self.tilezip), mode='r') as zf:
            with zf.open(self.tilezip.name[0:4].upper()+'.gml') as xmlf:
                txml=xxml.parse(xmlf).getroot()
        ginfo=_xmlgetchild(txml, ('member', 'ElevationGridCoverage', 'rectifiedGridDomain', 'RectifiedGrid'))
        return (int(_xmlgetchild(ginfo, 'offsetVector:1').text.split(' ')[0]), 
                int(_xmlgetchild(ginfo, 'offsetVector:0').text.split(' ')[1]))

    def opentile(self, skiplines):
        if self.zf is None and self.demf is None:
            self.zf=zipfile.ZipFile(str(self.tilezip), mode='r')
#            print('opened', str(self.tilezip))
            self.demf=self.zf.open(self.tilezip.name[0:4].upper()+'.asc')
            for _ in range(skiplines+5):
                self.demf.readline()
            return self.demf.readline
        else:
            raise RunTimeError('attempt to open elevations file in %s when already open' % self.tilezip)

    def closetile(self):
        if self.zf is None or self.demf is None:
            raise RunTimeError('attempt to close elevations file in %s when not open' % self.tilezip)
        else:
            self.demf.close()
            self.demf=None
            self.zf.close()
            self.zf=None
#            print('closed', str(self.tilezip))
