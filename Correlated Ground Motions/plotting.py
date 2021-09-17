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

    else:    
        intensity = 'None'
         
    palette = 'viridis'

##### Plot GMM Median Ground Motions #####    
    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log10(np.concatenate([np.linspace(.001,.01,9,endpoint=False),np.linspace(.01,.1,9,endpoint=False),np.linspace(.1,1.1,10,endpoint=False)]))
    cbar_label = ['.001']+8*['']+['.01']+8*['']+['.1']+8*['']+['1.0']      
    map = ax.imshow(np.log10(variables['data']/100),extent=extent,vmin=-3,vmax=0,origin='upper',cmap=palette)
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('Median GMM Estimated Ground Motion',pad = 30,fontsize=22)
    ch=plt.colorbar(map, shrink=0.7, ticks = cbar_vals)     
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('PGA (g)',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('sm_mean.png', dpi=300, bbox_inches='tight')
    plt.show(map)    
    


    
##### Plot Total Standard Deviation  ####   
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    map = ax.imshow(np.sqrt(variables['uncertaintydata'][0]**2 + variables['uncertaintydata'][1]**2) ,extent=extent,vmin=.25, vmax=.6, origin='upper',cmap=palette)   
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    th = plt.title('Total Residual Standard Deviation',loc='left',pad = 30,fontsize=16)
    ch=plt.colorbar(map, shrink=0.7)
    ch.ax.set_ylabel('Std of log PGA',fontsize=14)
    ch.ax.get_yaxis().labelpad = 2
    #plt.savefig('total_std.png', dpi=300, bbox_inches='tight')
    plt.show(map)    

##### Plot Within Standard Deviation  ####   
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    map = ax.imshow(variables['uncertaintydata'][0] ,extent=extent,vmin=.1, vmax=.6, origin='upper',cmap=palette)   
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    th = plt.title('Within-Event Residual Standard Deviation',loc='left',pad = 30,fontsize=16)
    ch=plt.colorbar(map, shrink=0.7)
    ch.ax.set_ylabel('Std of log PGA',fontsize=14)
    ch.ax.get_yaxis().labelpad = 2
    #plt.savefig('within_std.png', dpi=300, bbox_inches='tight')
    plt.show(map)
    

#### Plot Between Standard Deviation ####    
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    map = ax.imshow(variables['uncertaintydata'][1],extent=extent,vmin=.1, vmax=.6, origin='upper',cmap=palette) 
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    th = plt.title('Between-Event Residual Standard Deviation',loc='left',pad = 30,fontsize=16)
    ch=plt.colorbar(map, shrink=0.7)
    ch.ax.set_ylabel('Std of log PGA',fontsize=14)
    ch.ax.get_yaxis().labelpad = 2
    #plt.savefig('between_std.png', dpi=300, bbox_inches='tight')
    plt.show(map)


#### Plot Between Realization #####    
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log(np.array([1/8,1/4,1/2,1,2,4,8]))
    cbar_label = ['1/8','1/4','1/2','1','2','4','8']
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
    ch=plt.colorbar(map, shrink=0.7, ticks = cbar_vals)      
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('Multiplicative Factor',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('between_realization.png', dpi=300, bbox_inches='tight')
    plt.show(map)


#### Plot Within Realization ######
    fig = plt.figure(figsize=(10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log(np.array([1/8,1/4,1/2,1,2,4,8]))
    cbar_label = ['1/8','1/4','1/2','1','2','4','8']
    map = ax.imshow(out['cor']*variables['uncertaintydata'][0],extent=extent,vmin=np.log(1/8),vmax=np.log(8), origin='upper',cmap='seismic')
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'og',markersize=3,mew=.5,markerfacecolor = 'g',markeredgecolor='g',label='Reported Intensity')
        plt.plot(in_station_lons, in_station_lats, '^r',markersize=4,mew=.5,markerfacecolor = 'r',markeredgecolor='r',label='Seismic Instrument')
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('exp(Within Residual Realization)',pad = 30,fontsize=22)
    ch=plt.colorbar(map, shrink=0.7,ticks = cbar_vals)
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('Multiplicative Factor',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('within_realization.png', dpi=300, bbox_inches='tight')
    plt.show(map) 

        
#### Plot Ground Motion Realization #######        
    fig = plt.figure(figsize = (10,10))
    proj = cartopy.crs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.gridlines(crs=proj, draw_labels=True, linestyle='-')
    ax.set_extent(extent)
    cbar_vals = np.log10(np.concatenate([np.linspace(.001,.01,9,endpoint=False),np.linspace(.01,.1,9,endpoint=False),np.linspace(.1,1.1,10,endpoint=False)]))
    cbar_label = ['.001']+8*['']+['.01']+8*['']+['.1']+8*['']+['1.0']
    map = ax.imshow(np.log10(out['data_new']/100),extent=extent,vmin=-3,vmax=0, origin='upper',cmap=palette)
    if stationdata != 'None':
        plt.plot(sm_station_lons, sm_station_lats, 'g>', markersize = 4)
        plt.plot(in_station_lons, in_station_lats, 'r^', markersize = 6)
    ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'),facecolor=WATER_COLOR,alpha = 1,zorder=11)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linestyle=':',zorder=10)
    plt.plot(epicenter_lon,epicenter_lat,'*k',markersize = 15,label='Epicenter',alpha=.9,zorder = 12)
    plt.legend(loc='upper right',fontsize = 16)
    plt.title('Ground Motion Realization',pad = 30,fontsize=22)
    ch=plt.colorbar(map, shrink=0.7, ticks = cbar_vals)     
    ch.ax.set_yticklabels(cbar_label)
    ch.ax.set_ylabel('PGA (g)',fontsize=18)
    ch.ax.get_yaxis().labelpad = -1
    #plt.savefig('sm_realization.png', dpi=300, bbox_inches='tight')
    plt.show(map)


    return
