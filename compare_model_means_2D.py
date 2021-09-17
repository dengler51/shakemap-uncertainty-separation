#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:54:28 2020

@author: emthompson, dengler
"""

"""
Compare two the ratios of two shakemaps
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from mapio.shake import ShakeGrid
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import json
from impactutils.colors.cpalette import ColorPalette
from numpy import log10

SIZE = (5, 4)


def get_parser():
    desc = '''Compare two shakemaps.
This program is to quickly make maps comparing two different shakemaps. 
the arguments are paths to grid.xml files.
Note that ratios are grid1/grid2.
'''
    parser = argparse.ArgumentParser(description=desc, epilog='\n\n')
    parser.add_argument('grid1', type=str,
                        help='Path to a ShakeMap grid.xml file.')
    parser.add_argument('grid2', type=str,
                        help='Path to a ShakeMap grid.xml file.')
    parser.add_argument('-i', '--imt', type=str, default='pga',
                        help='Which IMT to use? A String such as pga, pgv, '
                             'psa03, psa10.')
    parser.add_argument('-l1', '--label1', type=str, default='Grid 1',
                    help='Label of quantity being plotted')
    parser.add_argument('-l2', '--label2', type=str, default='Grid 2',
                help='Label of quantity being plotted')
    parser.add_argument('-d','--data',type=str,default='None',
                        help='Folder with station location if desired')
    parser.add_argument('-o', '--output', type=str, default='compare.png',
                        help='Output filename. Default: compare.png')
    parser.add_argument('-n', '--nocoasts', action='store_true',
                        help='Suppresses printing of coastlines on the maps.')
    return parser


def main(pparser, args):
    """
    Args:
        pparser: ArgumentParser object.
        args: ArgumentParser namespace with command line arguments.
    """
    STATIONS_ZORDER = 1150
    
    g1 = ShakeGrid.load(args.grid1).getData()[args.imt]
    g2 = ShakeGrid.load(args.grid2).getData()[args.imt]
    
    header = ShakeGrid.load(args.grid1).getEventDict()
    
    fill = True
   
    ep_lat = header['lat']
    ep_lon = header['lon']
    
    label1 = args.label1
    label2 = args.label2
    
    if args.data != 'None':
        with open(args.data,'r') as f:
            stations = json.load(f)
        
        dimt = args.imt.lower()
        dimt = 'pga'
        
        # get the locations and values of the MMI observations
        mmi_dict = {
            'lat': [],
            'lon': [],
            'mmi': []
        }
        inst_dict = {
            'lat': [],
            'lon': [],
            'mmi': [],
            dimt: []
        }
        # Get the locations and values of the observed/instrumented
        # observations.
        for feature in stations['features']:
            lon, lat = feature['geometry']['coordinates']
            itype = feature['properties']['instrumentType'].lower()
            # If the network matches one of these then it is an MMI
            # observation
            if itype == 'observed':
                # Append data from MMI features
                amplitude = feature['properties']['intensity']
                if amplitude == 'null':
                    amplitude = 'nan'
                mmi_dict['mmi'].append(float(amplitude))
                mmi_dict['lat'].append(lat)
                mmi_dict['lon'].append(lon)
            else:
                # Otherwise, the feature is an instrument
    
                # If dimt is MMI then we have to use the converted value
                # from an instrumental value
                try:
                    mmi_conv = float(feature['properties']['intensity'])
                except ValueError:
                    mmi_conv = np.nan
                if dimt == 'mmi':
                    inst_dict[dimt].append(mmi_conv)
                    inst_dict['lat'].append(lat)
                    inst_dict['lon'].append(lon)
                elif len(feature['properties']['channels']) > 0:
                    # Other IMTs are given in the channel for instruments
                    channel = feature['properties']['channels'][0]
                    for amplitude in channel['amplitudes']:
                        if amplitude['name'] != dimt:
                            continue
                        if amplitude['value'] == 'null':
                            continue
                        inst_dict[dimt].append(float(amplitude['value']))
                        inst_dict['lat'].append(lat)
                        inst_dict['lon'].append(lon)
                        # Add mmi also for the fill color of symbols
                        inst_dict['mmi'].append(mmi_conv)
                        
                 
            
        

    g1_geodict = g1.getGeoDict()
    g2_geodict = g2.getGeoDict()
    try:
        cutdict = g1_geodict.getBoundsWithin(g2_geodict)
    except Exception:
        cutdict = g2_geodict.getBoundsWithin(g1_geodict)
    c1 = g1.interpolateToGrid(cutdict)
    c2 = g2.interpolateToGrid(cutdict)
    
    extent = (g1_geodict.xmin,g1_geodict.xmax,g1_geodict.ymin,g1_geodict.ymax)
    

    WATER_COLOR = [.47,.60,.81]
    
    a1 = c1.getData()
    a2 = c2.getData()
    ratio = np.log10(a1/a2)
   
    fig = plt.figure(figsize=SIZE)
    wid = 1.0
    height = 0.8
    
    vmin = -0.050
    vmax =  0.050

    fig = plt.figure(figsize=SIZE)
    
    proj = ccrs.PlateCarree()
    ax1 = plt.axes([x1, y1, wid, height], projection=proj)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    cs1 = ax1.imshow(np.flipud(ratio),extent=extent,vmin=vmin,vmax=vmax,
                      cmap=cmap)
    gl = ax1.gridlines(crs=proj, draw_labels=True, linestyle='-')
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'size': 8, 'color': 'black'}
    gl.top_labels = False
    gl.left_labels = False
    plt.plot(ep_lon,ep_lat,'*k',markersize = 10,label='Epicenter',alpha=.9)
    if args.data != 'None':
        if fill:
            intensity_colormap = ColorPalette.fromPreset('mmi')
            # Fill the symbols with the color for the intensity
            for i in range(len(mmi_dict['lat'])):
                mlat = mmi_dict['lat'][i]
                mlon = mmi_dict['lon'][i]
                mmi = mmi_dict['mmi'][i]
                mcolor = intensity_colormap.getDataColor(mmi)
                if np.isnan(mmi):
                    mcolor = (mcolor[0], mcolor[1], mcolor[2], 0.0)
                if i == 0:
                    ax1.plot(mlon, mlat, 'o', markerfacecolor='white',
                        markeredgecolor='k', markersize=2, mew=0.5,
                        zorder=STATIONS_ZORDER, transform=proj,label='Reported Intensity')
                else:    
                    ax1.plot(mlon, mlat, 'o', markerfacecolor='white',
                            markeredgecolor='k', markersize=2, mew=0.5,
                            zorder=STATIONS_ZORDER, transform=proj)
    
            for i in range(len(inst_dict['lat'])):
                mlat = inst_dict['lat'][i]
                mlon = inst_dict['lon'][i]
                #
                # TODO: Make the fill color correspond to the mmi
                # obtained from the IMT.
                #
                mmi = inst_dict['mmi'][i]
                mcolor = intensity_colormap.getDataColor(mmi)
                if np.isnan(mmi):
                    mcolor = (mcolor[0], mcolor[1], mcolor[2], 0.0)
                if i == 0:    
                    ax1.plot(mlon, mlat, '^',
                            markerfacecolor='white', markeredgecolor='k',
                            markersize=3, zorder=STATIONS_ZORDER, mew=0.5,
                            transform=proj,label='Seismic Instrument')
                else:
                    ax1.plot(mlon, mlat, '^',
                        markerfacecolor='white', markeredgecolor='k',
                        markersize=3, zorder=STATIONS_ZORDER, mew=0.5,
                        transform=proj)
        else:
            # Do not fill symbols
    
            # plot MMI as small circles
            mmilat = mmi_dict['lat']
            mmilon = mmi_dict['lon']
            ax1.plot(mmilon, mmilat, 'ko', fillstyle='none', mew=0.5,
                    markersize=4, zorder=STATIONS_ZORDER, transform=proj,label ='Reported Intensity')
    
            # plot instruments as slightly larger triangles
            instlat = inst_dict['lat']
            instlon = inst_dict['lon']
            ax1.plot(instlon, instlat, 'k^', fillstyle='none', mew=0.5,
                    markersize=6, zorder=STATIONS_ZORDER, transform=proj,label ='Seismic Instrument')
    
            plt.legend(loc ='upper right')
        
    
    if not args.nocoasts:
        ax1.add_feature(cfeature.COASTLINE.with_scale('50m'))
        ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle='-')
        ax1.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,zorder=11)
   
    plt.title('      Ratio of Proposed and W2018 Conditional Median PGA',pad=20,fontsize=14)
    ax_cbar1 = plt.axes([1.0, .2, 0.05, 0.8])
    
    cbar_vals = [log10(.9),log10(.92),log10(.94),log10(.96),log10(.98),log10(1),log10(1.02),log10(1.04),log10(1.06),log10(1.08),log10(1.1)]
    cbar_label = ['.9','.92','.94','.96','.98','1.0','1.02','1.04','1.06','1.08','1.1']
        
    cbar1 = fig.colorbar(cs1, cax=ax_cbar1,ticks=cbar_vals,
                        orientation='vertical')
    cbar1.ax.set_yticklabels(cbar_label)
    cbar1.ax.set_ylabel('Ratio (Proposed/W2018)',fontsize=12)
    cbar1.ax.get_yaxis().labelpad = 15

    # Difference plot
    ax1.set_xlim(cutdict.xmin, cutdict.xmax)
    ax1.set_ylim(cutdict.ymin, cutdict.ymax)

   plt.savefig(args.output, dpi=900, bbox_inches='tight')


if __name__ == '__main__':
    parser = get_parser()
    #pargs = parser.parse_args()
    
    pargs = parser.parse_args(['/Input Data/grid_alt.xml',
                                '/Input Data/grid_W2018.xml',
                                '-i', 'pga', '-l1', 'Proposed Conditioning Model Median PGA' ,'-l2', 'W2018 Conditioning Model Median PGA',
                                '-o','compare_models_mean.png' ,
                                '-d','/Input Data/stationlist.json'])
    main(parser, pargs)
