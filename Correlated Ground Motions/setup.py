import numpy as np
import math
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point
import time
from xml.dom import minidom

def initialize(shakemap, unc_INTRA, unc_INTER, stationdata, dm=1, dn=1,range_lat=[0,-1],range_lon=[0,-1]):
    """
    Set up grid spacing, site collections, and other data values
    
    INPUTS: 
    shakemap- shake grid of grid.xml
    uncertainty- shake grid of uncertainty.xml
    stationdata- data from stationlist.xml
    vscorr - boolean, determines if Vs30 are correlation. See
    JB2009
    dm, dn- vertical and horozontal spacing, respectively. Default value is 1
    OUT:
    Returns a dictionary with the following keys:    
    M,N - number of points vertically and horozontally
    K- number of stations
    uncertaintydata- array of uncertainty data at each point
    site_collection_SM- site collection for ShakeMap data
    site_collection_station- site collection for station data
    location_lat/lon_g- lat and lons for grid data, for plotting
    data- ShakeMap data at each point
    itensity- factor for non-native data
    site_collection_both- site collection for both SM and stations
    """

    start = time.time()

    attributes = shakemap.getAttributes()
    
    # Determine the size of the grid                                                                
    m = attributes['grid_specification']['nlat']
    n = attributes['grid_specification']['nlon']

    # M,N number of vertical, horizontal points considered
    if (range_lat[0]==0) & (range_lat[1]==-1):
        M = int(math.floor(m/dm))
    else:
        if range_lat[1]>0:
            M = int(math.floor((range_lat[1]-range_lat[0])/dm))
        else:
            M = int(math.floor((m-range_lat[1]+1-range_lat[0])/dm))

    if (range_lon[0]==0) & (range_lon[1]==-1):
        N = int(math.floor(n/dn))
    else:
        if range_lon[1]>0:
            N = int(math.floor((range_lon[1]-range_lon[0])/dn))
        else:
            N = int(math.floor((n-range_lon[1]+1-range_lon[0])/dn))


    # Determine number of stations
    if stationdata == 'None':
        K = 0
    else:
        K = np.size(stationdata['lat'])
    
    if unc_INTER == 'None':
        bias = True
    else:
        bias = False
#        bias = ComputeBias('Inputs/uncertainty2.xml')

    if bias == True:
        uncertainty = unc_INTRA.griddata
    else:
        rand = np.random.randn()
        inter_percent = .01
        
        #unc_INTER = inter_percent*unc_INTRA.griddata
        unc_INTER = unc_INTER.griddata
        unc_INTRA = unc_INTRA.griddata
#        if abs(rand) > 3:
#            rand = np.sign(rand)*(6 - np.sign(rand)*rand)
#        uncertainty = np.sqrt(np.power(rand*unc_INTER.griddata,2) + np.power(unc_INTRA.griddata,2))
        uncertainty = [unc_INTRA,unc_INTER]
        

    if (range_lon[0]!=0) | (range_lon[1]!=-1):
        if len(uncertainty)>2:
            uncertainty = uncertainty[range_lat[0]:range_lat[1],range_lon[0]:range_lon[1]]
        else:
            uncertainty[0] = uncertainty[0][range_lat[0]:range_lat[1],range_lon[0]:range_lon[1]]
            uncertainty[1] = uncertainty[1][range_lat[0]:range_lat[1],range_lon[0]:range_lon[1]]

    if len(uncertainty)>2:
        uncertainty = uncertainty[::dm,::dn]
    else:
        uncertainty[0] = uncertainty[0][::dm,::dn]
        uncertainty[1] = uncertainty[1][::dm,::dn]
        
#    uncertainty = uncertainty
#    uncertainty = uncertainty[:,:-1]
    # Allocate matrices for storing data and locations
    DATA = np.empty([M,N])
    location_lat_g = np.empty([M,N])
    location_lon_g = np.empty([M,N])
    site_SM = []

    # Puts griddata into Site class, then turns the sites into a SiteCollection
    for i in range(0,M):
        for j in range(0,N):

            DATA[i,j] = shakemap.griddata[range_lat[0]+i*dm,range_lon[0]+j*dn] # UNITS pctg  
            lat,lon = shakemap.getLatLon(range_lat[0]+i*dm,range_lon[0]+j*dn)
            location_lat_g[i,j] = lat
            location_lon_g[i,j] = lon
            site_SM.append(Site(location = Point(location_lon_g[i,j], location_lat_g[i,j]), 
                                    vs30 = 760, vs30measured = True, z1pt0 = 100, z2pt5 = 1))
    
    site_collection_SM = SiteCollection(site_SM)

    # Store lat lon points
    location_lat_s = np.empty([K])
    location_lon_s = np.empty([K])
    intensity = np.empty([K])
    site_station = []

    # Puts stationdata into Site class, then turns the sites into a SiteCollection
    if stationdata != 'None':
        for i in range(0,K):
            location_lat_s[i] = stationdata['lat'][i]
            location_lon_s[i] = stationdata['lon'][i]
            if stationdata['name'][i] == 'DERIVED':
                intensity[i] = 1
            else:
                intensity[i] = 0
        
            site = Site(location=Point(location_lon_s[i], location_lat_s[i]),
                        vs30 = 760, vs30measured = True, z1pt0 = 100, z2pt5 = 1)
            site_station.append(site)

        site_collection_station = SiteCollection(site_station)
        site_both = site_station + site_SM
    else:
        site_both = site_SM
        site_collection_station = 'None'
        intensity = 'None'
        
    site_collection_both = SiteCollection(site_both)


    
    end = time.time() - start
    print('%6f Initialization Time:' %(end))
    
    return {'M': M, 'N': N, 'K': K, 'uncertaintydata':uncertainty, \
                'site_collection_SM':site_collection_SM, \
                'site_collection_station':site_collection_station, \
                'location_lat_g': location_lat_g, 'location_lon_g':location_lon_g, \
                'data':DATA, 'intensity':intensity, \
                'site_collection_both':site_collection_both}

def ComputeBias(uncertainty):
    root = minidom.parse(uncertainty)
    evu = root.getElementsByTagName('event_specific_uncertainty')
    for instance in evu:
        value = instance.getAttribute('value')
        if float(value) == -1.0:
            return False
        else:
            return True

            
