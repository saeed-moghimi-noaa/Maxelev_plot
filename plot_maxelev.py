#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot maxele

Maxeles for ops, para:

    *GESTOFS operational 20210406t00z: https://noaa-gestofs-pds.s3.amazonaws.com/estofs.20210406/estofs.t00z.fields.cwl.maxele.nc
    *GESTOFS parallel 20210406t00z: http://polar.ncep.noaa.gov/estofs/glo_para/fields/estofs.20210406/estofs.t00z.fields.cwl.maxele.nc

Other maxeles that could be useful to plot:

     Pre-Boston crash/cold start fix: https://polar.ncep.noaa.gov/estofs/autoval/tests/estofs.20210202.crash/estofs.20210202/estofs.t00z.fields.cwl.maxele.nc   
    Guoming's 2018010-201812 GESTOFS v4 grid: https://drive.google.com/file/d/1lJw6P8zpJ9QpPSf7W1gdeGg7khTVmbv2/view?usp=drive_web

List of spots: https://docs.google.com/spreadsheets/d/16Hw0bwdBytrPloIu_oJs11VBBUGlw9QqLNVlTojJWaA/edit#gid=0


"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2021, UCAR/NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"



from   matplotlib.colors import LinearSegmentedColormap
import netCDF4
import cartopy.crs as ccrs
import matplotlib.tri as Tri
import matplotlib.pyplot as plt
import sys,os
from geo_regions import get_region_extent
from geo_regions import defs
import numpy as np
import shapefile   
from glob import glob
from math import floor
from matplotlib import patheffects
from matplotlib.image import imread

#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.io.img_tiles as cimgt
#from cartopy.io.img_tiles import GoogleTiles

# class ShadedReliefESRI(GoogleTiles):
#     # shaded relief
#     def _image_url(self, tile):
#         x, y, z = tile
#         url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
#                'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
#                z=z, y=y, x=x)
#         return url

def ReadTri(fmaxele):

    """
    fname: one of fort.*.nc file
    tind: time index
    """ 
    fname =  os.path.abspath(fmaxele)
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    x   = nc.variables['x'][:]
    y   = nc.variables['y'][:]
    d   = nc.variables['depth'][:]
    zeta = nc.variables['zeta_max'][:].squeeze()
    #zeta = np.ma.masked_where(np.isnan(zeta),zeta)
    # read connectivity array
    el  = nc.variables['element'][:] - 1
    # create a triangulation object, specifying the triangle connectivity array
    print ('[info:] Generate Tri ...')
    tri  = Tri.Triangulation(x,y, triangles=el)
    nc.close()
    return x,y,tri,d,zeta

def make_map(projection=ccrs.PlateCarree(), xylabels = False):
    
    """
    Generate fig and ax using cartopy
    input: projection
    output: fig and ax
    """

    subplot_kw = dict(projection=projection)
    fig, ax = plt.subplots(figsize=(9, 13),
                           subplot_kw=subplot_kw)
    if xylabels:
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    else:
        gl = ax.gridlines(draw_labels=False)
        gl.xlines = False
        gl.ylines = False

    return fig, ax

def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return floor( ( lon + 180 ) / 6) + 1


def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000):
    """

    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    # Projection in metres
    utm = ccrs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, zorder=3)


############################
cdict = {'red': ((0.  , 1  , 1),
                 (0.05, 1  , 1),
                 (0.11, 0  , 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1   , 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

jetMinWi = LinearSegmentedColormap('my_colormap',cdict,256)
cmap = jetMinWi

data_dir = '/mnt/c/Users/Saeed.Moghimi/Documents/work/linux_working/02-models/02-adcirc/02-maxelev/data/'


fv1_opr           = data_dir + 'estofs.t00z.fields.cwl.maxele_current_operational.nc'
title_opr         = 'GESTOFS-Operational'
fv4_para          = data_dir + 'estofs.t00z.fields.cwl.maxele_v4_upgrade.nc'
title_v4          = 'GESTOFS-Upgrade'

fund_3month_barry = data_dir + 'maxele_UND_3month_GESTOFS_v4.63.nc'
title_v4_und          = 'GESTOFS-Upgrade-barry'


fv1_opr_c           = data_dir + 'tampa_bay_high_water_estofs.t00z.fields.cwl.maxele.nc'
title_v1_c          = 'GESTOFS-Operational-C'

fnames = [fv1_opr,   fv4_para, fv1_opr_c ]#, fund_3month_barry]
titles = [title_opr, title_v4, title_v1_c ]#, title_v4_und]

###
shpFilePath = '/home/moghimis/linux_working/02-models/02-adcirc/01-meshes/04-GESTOFS-shapefiles/'  
shplist = glob( shpFilePath + '*.shp')

#regions = ['ike_local', 'ike_region', 'san_delaware','san_jamaica_bay','Tampa_Area_m', 'Marshall','Palau','and_local_lu']    
#regions = ['san_newyork','Port_Arthur_m',  'Tampa_Area_m','ike_local','ike_region']
#regions = ['Palau','and_local_lu'] 
regions = ['Tampa_Area_m','san_jamaica_bay'] 




bkg_img = '/home/moghimis/linux_working/opt/NE1_50M_SR_W/NE1_50M_SR_W.tif'

#regions = []
#regions.append( 'hsofs_region')

#regions.append( 'and_local_lu')
#regions = ['puertorico','isa_local','san_newyork','san_delaware','san_jamaica_bay',
#    'and_fl_lu','NYC_area','Tampa_area','Marshall','Palau',
#    'ike_region','ike_local','ike_wave']
#regions.append( 'hsofs_region')
#regions.append( 'and_local_lu')
#regions.append( 'san_area2')
#regions.append( 'san_track')


for ifname in range(len(fnames)):
    fname = fnames[ifname]
    print ('> Maxelev file: ', fname)
    for region in regions:
        print ('  > Maxelev region: ', region)

        defs['elev']['vmin']  = 0.2
        defs['elev']['vmax']  = 2.2

        ### construct levels
        var    = defs['elev']
        vmin   = var['vmin']
        vmax   = var['vmax']
        dv     = (vmax-vmin)/21.0
        levels = np.arange(vmin,vmax+dv,dv)

        ### read data #########
        lon,lat,tri,dep,val  = ReadTri(fname)

        #stamen_terrain = cimgt.Stamen('terrain-background')
        #fig, ax = make_map( projection=stamen_terrain.crs)
        #fig, ax = make_map(projection=ShadedReliefESRI().crs)
        #ax.add_image(ShadedReliefESRI(), 10, zorder=-1)

        fig, ax = make_map()
        fig.set_size_inches(9,9)

        lim = get_region_extent(region = region)
        extent = [lim['xmin'],lim['xmax'],lim['ymin'],lim['ymax']]       
        ax.set_extent(extent)
        ax.set_title(titles[ifname])
        
        source_proj = ccrs.PlateCarree()
        ax.imshow(imread(bkg_img), origin='upper', transform=source_proj, 
          extent=extent,zorder=0)

        #sys.exit()

        ### maxelev
        if True:
            cf1    = ax.tricontourf(tri,val.data,levels=levels, cmap = cmap , extend='both',alpha = 1.0,zorder=2)
            cb     = plt.colorbar(cf1,shrink = 0.3,ticks = [vmin,(vmin+vmax)/2,vmax])   
            cb.set_label('TWL [m]')
        else:
            ax.tripcolor(tri, val.data, cmap = cmap, clim = (vmin,vmax))
            plt.colorbar()
        ### plot mesh
        if False:
            ax.triplot(tri, 'k-', lw=0.1, alpha=0.5)

        
        ### shoreline
        if True:
            if False:
                cond1 = ax.tricontour(tri,dep ,levels=[0.0]  ,colors='k',linewidths=[0.5], alpha= 1.0,zorder=3)
            else:
                print ('[info:] Read shapefiles ...')
                for fshp in shplist[:]:
                    sf = shapefile.Reader(fshp)
                    for shape in sf.shapeRecords():
                        x = [i[0] for i in shape.shape.points[:]]
                        y = [i[1] for i in shape.shape.points[:]]
                        ax.plot(x,y,'k',lw=0.2)

            
        if True:
            scale_bar(ax, ax.projection, 10) 
            
        filename = 'figs3/maxelev_' + region + '_'+ titles[ifname]
        filename = filename.replace('.','-')
        plt.savefig(filename + '.png',dpi=600)
        plt.close('all')

