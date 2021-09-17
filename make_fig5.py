#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 10:54:28 2020

@author: emthompson
"""

"""
Compare two shakemaps.
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

SIZE = (5, 4)


def get_parser():
    desc = '''Compare two shakemaps.
This program is to quickly make maps comparing two different shakemaps. Since
the main goal is to compare across maps using ShakeMap 3.5 and ShakeMap 4.0
the arguments are paths to grid.xml files.
Note that ratios are grid1/grid2 and differences are grid1 - grid2.
'''
    parser = argparse.ArgumentParser(description=desc, epilog='\n\n')
    parser.add_argument('grid', type=str,
                        help='Path to a ShakeMap grid.xml file.')
    parser.add_argument('-i', '--imt', type=str, default='pga',
                        help='Which IMT to use? A String such as pga, pgv, '
                             'psa03, psa10.')
    parser.add_argument('-l', '--label', type=str, default='Grid 1',
                    help='Label of quantity being plotted')
    parser.add_argument('-s','--stat',type=str,default='mean')
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
    
    g = ShakeGrid.load(args.grid).getData()[args.imt]
    
    header = ShakeGrid.load(args.grid).getEventDict()
    
    fill = True
   
    ep_lat = header['lat']
    ep_lon = header['lon']
    
    label = args.label
       
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
                        
                 
            
        

    g_geodict = g.getGeoDict()
    
    c = g.interpolateToGrid(g_geodict)
   
    WATER_COLOR = [.47,.60,.81]
    
    a = c.getData()
    
    if (args.imt.lower() == 'pga') & ((args.stat == 'mean') | (args.stat =='median') ):
        a = np.log10(a/100)
 
    lats = np.linspace(g_geodict.ymin, g_geodict.ymax, a.shape[0])
    lons = np.linspace(g_geodict.xmin, g_geodict.xmax, a.shape[1])

    fig = plt.figure(figsize=SIZE)
    wid = 1.0
    height = 0.8

    # Ratio plot
    # color_range = eval(args.range)
    # levels = list(np.linspace(color_range[0], color_range[1], int(color_range[1]-color_range[0]+1)))
    # levels = list(np.linspace(color_range[0], color_range[1], 256))
    #levels = list(np.linspace(np.log10(.92), np.log10(1.08),12))
    #cmap = plt.cm.Spectral_r
    
    cmap = 'viridis'
    
    extent = (g_geodict.xmin,g_geodict.xmax,g_geodict.ymin,g_geodict.ymax)

    x1 = 0.05
    y1 = 0.2
    
    if (args.imt=='pga') & (args.stat == 'mean'):
        vmin = -3
        vmax = np.log10(1)


    else:
        vmin=0
        vmax=.2
        vmin=.25
        vmax=.6

    fig = plt.figure(figsize=SIZE)
    
    proj = ccrs.PlateCarree()
    ax1 = plt.axes([x1, y1, wid, height], projection=proj)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    cs1 = ax1.imshow(np.flipud(a),extent=extent,vmin = vmin, vmax = vmax,
                      cmap=cmap)
    # cs1 = ax1.contourf(lons, lats, np.flipud(np.log10(diff2)),
    #                    cmap=cmap, extend='both')
    gl = ax1.gridlines(crs=proj, draw_labels=True, linestyle='-')

    #gl.xformatter = LongitudeFormatter()
    #gl.yformatter = LatitudeFormatter()
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
                    ax1.plot(mlon, mlat, 'o', markerfacecolor=mcolor,
                        markeredgecolor='k', markersize=4, mew=0.5,
                        zorder=STATIONS_ZORDER, transform=proj,label='Reported Intensity')
                else:    
                    ax1.plot(mlon, mlat, 'o', markerfacecolor=mcolor,
                            markeredgecolor='k', markersize=4, mew=0.5,
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
                            markerfacecolor=mcolor, markeredgecolor='k',
                            markersize=6, zorder=STATIONS_ZORDER, mew=0.5,
                            transform=proj,label='Seismic Instrument')
                else:
                    ax1.plot(mlon, mlat, '^',
                        markerfacecolor=mcolor, markeredgecolor='k',
                        markersize=6, zorder=STATIONS_ZORDER, mew=0.5,
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
    
            # mmi_file = open(args.data+'/mmi_sta_locs.txt','r')
            # sta_file = open(args.data+'/inst_sta_locs.txt','r')
            # read_mmi_flag = 1
            # mmi_lat = []
            # mmi_lon = []
            # while read_mmi_flag == 1:
            #     line = mmi_file.readline()
            #     if line != '':
            #         line = eval('['+line.replace('\n','')+']')
            #         mmi_lat.append(line[0])
            #         mmi_lon.append(line[1])
            #     else:
            #         read_mmi_flag = 0
            # read_sta_flag = 1        
            # sta_lat = []
            # sta_lon = []
            # while read_sta_flag == 1:
            #     line = sta_file.readline()
            #     if line != '':
            #         line = eval('['+line.replace('\n','')+']')
            #         sta_lat.append(line[0])
            #         sta_lon.append(line[1])
            #     else:
            #         read_sta_flag = 0
        # plt.plot(mmi_lon,mmi_lat,'ok',markersize=2,mew=.5,markerfacecolor = 'w',markeredgecolor='k',label='Reported Intensity')
        # plt.plot(sta_lon,sta_lat,'^k',markersize=3,mew=.5,markerfacecolor = 'w',markeredgecolor='k',label='Seismic Instrument')
        # plt.xlim([np.min(lons),np.max(lons)])
        # plt.ylim([np.min(lats),np.max(lats)])
        plt.legend(loc='upper right')
        
    
    if not args.nocoasts:
        ax1.add_feature(cfeature.COASTLINE.with_scale('50m'))
        ax1.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle='-')
        ax1.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,zorder=11)
    if args.imt == 'pga':
        imt_str = r'$log_{10}$ PGA'
    elif args.imt == 'stdpga':
        imt_str = 'Std of PGA'
    #print(np.max(a))
    #print(np.min(a))
    plt.title(label,loc='left',pad = 30,fontsize=16)
    ax_cbar1 = plt.axes([1.0, .2, 0.05, 0.8])
    
    if args.stat == 'mean': 
        cbar_vals = np.log10(np.concatenate([np.linspace(.001,.01,9,endpoint=False),np.linspace(.01,.1,9,endpoint=False),np.linspace(.1,1.1,10,endpoint=False)]))
        #cbar_vals = [-2,np.log10(.02),np.log10(.04),np.log10(.08),-1,np.log10(.2),np.log10(.4),np.log10(.8),0]
        #cbar_label = ['.01','.02','.04','.08','.1','.2','.4','.8','1.0']
        cbar_label = ['.001']+8*['']+['.01']+8*['']+['.1']+8*['']+['1.0']
        
        cbar1 = fig.colorbar(cs1, cax=ax_cbar1,ticks = cbar_vals,
                        orientation='vertical')
        
        cbar1.ax.set_yticklabels(cbar_label)
    else:
        cbar1 = fig.colorbar(cs1, cax=ax_cbar1,
                        orientation='vertical')

    # cbar1 = fig.colorbar(cs1, cax=ax_cbar1,
    #                      orientation='horizontal',
    #                      ticks=levels)
   
    # cbar_label = []
    # for i in range(len(levels)):
    #     cbar_label.append(f'{levels[i]:.1f}')
    # #print(cbar_label)  
    # cbar_label = [r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$0.01$',r'$0.1$',r'$1$',r'$10$']
    # cbar1.ax.set_xticklabels(cbar_label)
    
    
    for t in cbar1.ax.get_xticklabels():
        t.set_fontsize(8)
#    cbar1.ax.set_xlabel('Point-Source MMI Uncertainty/\n Atlas MMI Uncertainty')
#    cbar1.ax.set_xlabel('Atlas  MMI /\n ShakeMap 3.5 MMI')
    #cbar1.set_label(r'$log_{10}$ PGA')
    if (args.imt == 'pga') & (args.stat == 'mean'):
        cbar1.ax.set_ylabel('PGA (g)',fontsize=14)
    else:
        cbar1.ax.set_ylabel('Std of log PGA',fontsize=14)
    cbar1.ax.get_yaxis().labelpad = 2

    # Difference plot
    levels = list(np.linspace(-.55, .55, 8))
    levels = list(np.linspace(-.1, .1, 12))
    x1 = 0.55
    ax1.set_xlim(g_geodict.xmin, g_geodict.xmax)
    ax1.set_ylim(g_geodict.ymin, g_geodict.ymax)

#    ax2 = plt.axes([x1, y1, wid, height], projection=ccrs.PlateCarree())
#    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
#    cs2 = ax2.contourf(lons, lats, np.flipud(dif), levels,
#                       cmap=cmap, extend='both')
#    if not args.nocoasts:
#        ax2.add_feature(cfeature.COASTLINE.with_scale('50m'))
#        ax2.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle='-')
#        ax2.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,zorder=11)
#
#    ax_cbar2 = plt.axes([x1, y1-0.1, wid, 0.05])
#    cbar2 = fig.colorbar(cs2, cax=ax_cbar2,
#                         orientation='horizontal',
#                         ticks=levels)
#    for t in cbar2.ax.get_xticklabels():
#        t.set_fontsize(4)
#    if args.imt == 'pgv':
#        cbar2.ax.set_xlabel('%s Difference (cm/s)' % args.imt)
#    elif args.imt == 'stdmmi':
##        cbar2.ax.set_xlabel('Point-Source MMI Uncertainty -\n Atlas MMI Uncertainty')
#        cbar2.ax.set_xlabel('Atlas (Original) MMI Uncertainty -\n Atlas (Fixed) MMI Uncertainty')
#    else:
#        cbar2.ax.set_xlabel('%s Difference (percent g)' % args.imt)
#    cbar2.ax.get_yaxis().labelpad = 15
    plt.savefig(args.output, format='jpg', dpi=1200, bbox_inches='tight')


if __name__ == '__main__':
    parser = get_parser()
    #pargs = parser.parse_args()
    
    pargs = parser.parse_args(['/Users/dengler/shakemap_profiles/default/data/ci38457511/Current/uncertainty_alt.xml',
                                '-i', 'stdpga','-s','std', '-l', 'Proposed Model Conditional Ground Motion', 
                                '-o', '/Users/dengler/Documents/Codes/Loss with Uncertainty/W2018_shakemap.jpeg',
                                ])    
    main(parser, pargs)