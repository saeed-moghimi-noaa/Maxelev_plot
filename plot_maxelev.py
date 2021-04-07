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

fnames = [fv1_opr,   fv4_para, fund_3month_barry]
titles = [title_opr, title_v4, title_v4_und]

###

#regions = ['ike_local', 'ike_region', 'san_delaware','san_jamaica_bay','Tampa_Area_m', 'Marshall','Palau','and_local_lu']    
regions = ['NYC_Area_m','Port_Arthur_m', 'san_newyork']


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

        fig, ax = make_map()
        fig.set_size_inches(9,9)

        lim = get_region_extent(region = region)
        extent = [lim['xmin'],lim['xmax'],lim['ymin'],lim['ymax']]       
        ax.set_extent(extent)
        ax.set_title(titles[ifname])

        ### shoreline
        cond1 = ax.tricontour(tri,dep ,levels=[0.0]  ,colors='k',linewidth=0.01, alpha= 0.7)
        ### maxelev
        if True:
            cf1    = ax.tricontourf(tri,val.data,levels=levels, cmap = cmap , extend='both',alpha = 1.0)#extend='max' )  #,extend='both'
            cb     = plt.colorbar(cf1,shrink = 0.3,ticks = [vmin,(vmin+vmax)/2,vmax])   
            cb.set_label('TWL [m]')
        else:
            ax.tripcolor(tri, val.data, cmap = cmap, clim = (vmin,vmax))
            plt.colorbar()
        ### plot mesh
        if False:
            ax.triplot(tri, 'k-', lw=0.1, alpha=0.5)
                            
                                                
        filename = 'figs/maxelev_' + region + '_'+ titles[ifname]
        filename = filename.replace('.','-')
        plt.savefig(filename + '.png',dpi=450)
        plt.close('all')


