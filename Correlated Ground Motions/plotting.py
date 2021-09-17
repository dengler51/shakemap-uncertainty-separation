import cartopy
import matplotlib.pyplot as plt
import numpy as np
import scipy
from neicio.gmt import GMTGrid
from matplotlib import cm
from scipy.stats import norm

WATER_COLOR = [.47,.60,.81]


def plot(out, variables, voi, shakemap, stationdata):
#    maxdata = np.amax(out['data_new'])
    extent = (np.min(variables['location_lon_g']),np.max(variables['location_lon_g']),np.min(variables['location_lat_g']),np.max(variables['location_lat_g']))
    
    attributes = shakemap.getAttributes()
    epicenter_lat = attributes['event']['lat']
    epicenter_lon = attributes['event']['lon']
    
    if len(stationdata) == 2:
        read_mmi_flag = 1
        in_station_lats = []
        in_station_lons = []
        while read_mmi_flag == 1:
            line = stationdata[0].readline()
            if line != '':
                line = eval('['+line.replace('\n','')+']')
                in_station_lats.append(line[0])
                in_station_lons.append(line[1])
            else:
                read_mmi_flag = 0
        read_sta_flag = 1        
        sm_station_lats = []
        sm_station_lons = []
        while read_sta_flag == 1:
            line = stationdata[1].readline()
            if line != '':
                line = eval('['+line.replace('\n','')+']')
                sm_station_lats.append(line[0])
                sm_station_lons.append(line[1])
            else:
                read_sta_flag = 0
#    elif stationdata != 'None':
#        intensity = stationdata['name'] 
#
#        SM = []
#        IN = []
#        for value in enumerate(intensity):
#            if ((value == 'UNCERTAINTY')or(value == 'DYFI')or(value == 'MMI')or(value == 'CIIM')):
#    
#                IN.append(value[0])
#            else:
#                SM.append(value[0])
#    
#        sm_station_lons = [stationdata['lon'][j] for j in SM]
#        sm_station_lats = [stationdata['lat'][j] for j in SM]
#        in_station_lats = [stationdata['lat'][j] for j in IN]
#        in_station_lons = [stationdata['lon'][j] for j in IN]  
#
#    
    else:    
        intensity = 'None'
         
    palette = 'viridis'
    #palette = cm.jet

    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
#    map = ax.imshow(out['cor']*variables['uncertaintydata'],extent=shakemap.getRange(),vmin=-3,vmax=3, origin='upper',cmap=palette)
    cbar_vals = np.log(np.array([1/8,1/4,1/2,1,2,4,8]))
    cbar_label = ['1/8','1/4','1/2','1','2','4','8']
    map = ax.imshow(out['cor']*variables['uncertaintydata'][0],extent=extent,vmin=np.log(1/8),vmax=np.log(8), origin='upper',cmap='seismic')
    #map = ax.imshow((norm.rvs(size=1)*variables['uncertaintydata'][0]),extent=extent,vmin=-3,vmax=3, origin='upper',cmap=palette)

    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'og',markersize=3,mew=.5,markerfacecolor = 'g',markeredgecolor='g',label='Reported Intensity')
        plt.plot(in_station_lons, in_station_lats, '^r',markersize=4,mew=.5,markerfacecolor = 'r',markeredgecolor='r',label='Seismic Instrument')
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('exp(Within Residual Realization)',pad = 30,fontsize=22)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('Within Event Residual Realization ('+voi+')\n for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7,ticks = cbar_vals)
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('Multiplicative Factor',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('within_realization.png', dpi=300, bbox_inches='tight')
    plt.show(map)
    
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
#    map = ax.imshow(out['cor']*variables['uncertaintydata'],extent=shakemap.getRange(),vmin=-3,vmax=3, origin='upper',cmap=palette)
    
    #map = ax.imshow(np.sqrt(variables['uncertaintydata'][0]**2 + variables['uncertaintydata'][1]**2),extent=extent,vmin=.0, vmax=.7, origin='upper',cmap=palette)
    map = ax.imshow(variables['uncertaintydata'][0] ,extent=extent,vmin=.2, vmax=.6, origin='upper',cmap=palette)
   
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)

    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    th = plt.title('Current Model Updated Total STD ('+voi+')\n for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.savefig('W2018_std.png', dpi=300, bbox_inches='tight')
    plt.show(map)
    
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
#    map = ax.imshow(out['cor']*variables['uncertaintydata'],extent=shakemap.getRange(),vmin=-3,vmax=3, origin='upper',cmap=palette)
    if voi == 'MMI':
       map = ax.imshow(variables['uncertaintydata'][1],extent=extent,vmin=0, vmax=.7, origin='upper',cmap=palette)
    else:
       map = ax.imshow(variables['uncertaintydata'][1],extent=extent,vmin=.2, vmax=.6, origin='upper',cmap=palette)
 
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)

    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('Updated Between-Event STD ('+voi+')\n for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.savefig('proposed_std.png', dpi=300, bbox_inches='tight')
    plt.show(map)

    
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log(np.array([1/8,1/4,1/2,1,2,4,8]))
    cbar_label = ['1/8','1/4','1/2','1','2','4','8']
#    map = ax.imshow(out['cor']*variables['uncertaintydata'],extent=shakemap.getRange(),vmin=-3,vmax=3, origin='upper',cmap=palette)
    if voi == 'MMI':
        map = ax.imshow(1.5*variables['uncertaintydata'][1],extent=extent,vmin=-2, vmax=2, origin='upper',cmap='seismic')
    else:
        map = ax.imshow(1.5*variables['uncertaintydata'][1],extent=extent,vmin=np.log(1/8),vmax=np.log(8), origin='upper',cmap='seismic')

    
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('exp(Between Residual Realization)',pad = 30,fontsize=22)
    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('Between Event Residual Realization ('+voi+')\n for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7, ticks = cbar_vals)
      
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('Multiplicative Factor',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('between_realization.png', dpi=300, bbox_inches='tight')
    plt.show(map)

    
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    map = ax.imshow(variables['uncertaintydata'][0],extent=extent,vmin=.1,vmax=.6, origin='upper',cmap=palette)
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 6)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)

    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('Shaking Within Uncertainty ('+voi+') for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    plt.show(map)

    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    map = ax.imshow(variables['uncertaintydata'][1],extent=extent,vmin=.1,vmax=.6, origin='upper',cmap=palette)
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'og',markersize=3,alpha=.4,mew=.5,markerfacecolor = 'g',markeredgecolor='g',label='Reported Intensity')
        plt.plot(in_station_lons, in_station_lats, '^r',markersize=3,alpha=.4,mew=.5,markerfacecolor = 'r',markeredgecolor='r',label='Seismic Instrument')
        plt.legend()
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)


    locstr = attributes['event']['event_description']
    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('Shaking Between Uncertainty ('+voi+') for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7)
    
    plt.show(map)


    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log10(np.concatenate([np.linspace(.001,.01,9,endpoint=False),np.linspace(.01,.1,9,endpoint=False),np.linspace(.1,1.1,10,endpoint=False)]))
        #cbar_vals = [-2,np.log10(.02),np.log10(.04),np.log10(.08),-1,np.log10(.2),np.log10(.4),np.log10(.8),0]
        #cbar_label = ['.01','.02','.04','.08','.1','.2','.4','.8','1.0']
    cbar_label = ['.001']+8*['']+['.01']+8*['']+['.1']+8*['']+['1.0']

    if voi == 'MMI':
        map = ax.imshow(variables['data'],extent=extent,vmin=-1, vmax=9.8, origin='upper',cmap=palette)
    else:
        map = ax.imshow(np.log10(variables['data']/100),extent=extent,vmin=-3,vmax=0,origin='upper',cmap=palette)

    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('Median GMM Estimated Ground Motion',pad = 30,fontsize=22)

    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('ShakeMap ('+voi+') for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7, ticks = cbar_vals)
      
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('PGA (g)',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    plt.savefig('sm_mean.png', dpi=300, bbox_inches='tight')
    plt.show(map)
    
    if voi == 'MMI':
        out['data_new'][out['data_new']<0]=0
        out['data_new'][out['data_new']>9]=9
        
    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log10(np.concatenate([np.linspace(.001,.01,9,endpoint=False),np.linspace(.01,.1,9,endpoint=False),np.linspace(.1,1.1,10,endpoint=False)]))
        #cbar_vals = [-2,np.log10(.02),np.log10(.04),np.log10(.08),-1,np.log10(.2),np.log10(.4),np.log10(.8),0]
        #cbar_label = ['.01','.02','.04','.08','.1','.2','.4','.8','1.0']
    cbar_label = ['.001']+8*['']+['.01']+8*['']+['.1']+8*['']+['1.0']
    
    if voi == 'MMI':
        map = ax.imshow(out['data_new'],extent=extent,vmin=-1, vmax=9.8, origin='upper',cmap=palette)
    else:
        map = ax.imshow(np.log10(out['data_new']/100),extent=extent,vmin=-3,vmax=0, origin='upper',cmap=palette)
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    locstr = attributes['event']['event_description']
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('Ground Motion Realization',pad = 30,fontsize=22)

    mag = attributes['event']['magnitude']
    datestr = attributes['event']['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
    #th = plt.title('Shakemap Realization ('+voi+') for %s - %s M%.1f, (epsilon)' % (locstr,datestr,mag), y = 1.08)
    ch=plt.colorbar(map, shrink=0.7, ticks = cbar_vals)
      
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('PGA (g)',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('sm_realization.png', dpi=300, bbox_inches='tight')
    plt.show(map)


    return
