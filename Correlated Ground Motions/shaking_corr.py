#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:32:38 2020

@author: sverros, dengler
"""

#######################################################
# Code for computing the spatial correlation for a ShakeMap,
# adding to a ShakeMap grid, and computing multiple realizations
# VARIABLES:
#     voi - variable of interest, i.e. PGA
#     r - radius of influence
#     num_realization- integer for desired number of realizations
#     corr_model- JB2009 or GA2010
#     vs_corr- Vs30 correlated bool, see JB2009
#     input data- grid.xml, uncertainty.xml, and stationlist.xml
#         stored in Inputs directory
#######################################################

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import truncnorm
from neicio.readstation import readStation
from neicio.shake import ShakeGrid
from neicio.gmt import GMTGrid
import time
from matplotlib import cm
import sys
from setup import initialize
from loop import main
from realizations import realizations
from plotting import plot

                                                                                                                                                                                                 
voi = 'PGA'
r = [15]
num_realizations = 0
cor_model = 'LB2013'
vscorr = True
plot_on = True

np.random.seed(88532) 
 

intra_flag = 1
station_flag = 0


mmi_file = '/../Input Data/grid_alt.xml'
phi_file = '/../Input Data/uncertainty_alt_within.xml'
tau_file = '/../Input Data/uncertainty_alt_between.xml'
#Default Range
range_lat = [0,-1]
range_lon = [0,-1]
##

for R in range(0,np.size(r)):
    radius = r[R]
    # Get shakemap for desired variable, PGA, uncertainty grid and stationdata
    shakemap = ShakeGrid(mmi_file, variable = '%s' % voi)

    # Uncertainty Data: Units in ln(pctg)
    if intra_flag == 1:
        unc_INTRA = ShakeGrid(phi_file, variable= 'STD%s' % voi)
        unc_INTER = ShakeGrid(tau_file, variable= 'STD%s' % voi)
    else:
        unc_INTRA = ShakeGrid(unc_file, variable= 'STD%s' % voi)
        unc_INTER = 'None'
        
    if station_flag == 1:
        # Station Data: Units in pctg
        stationlist = station_file
        stationdata = readStation(stationlist)
    else:
        stationlinst = 'None'
        stationdata = 'None'
        
        
    print('Calling initialize')
    variables = initialize(shakemap, unc_INTRA, unc_INTER,stationdata,dm=1,dn=1, range_lat=range_lat,range_lon=range_lon)

    print('Radius: ', radius)
    print(variables['K'], ' stations', variables['N']*variables['M'], ' points')

    rand = []    
    rand = np.random.randn(variables['N']*variables['M'])
    
    # Truncate to +/- 3 sigma
    for i in range(0,variables['N']*variables['M']):
        while abs(rand[i])>3:
            rand[i]=np.random.randn(1)

    print('Calling main')
    out = main(variables, r, voi, rand, cor_model, vscorr)
    
    if num_realizations > 0:
        print('Computing realizations')
        realizations(num_realizations, radius, variables['N'], variables['M'], out['grid_arr'],
                 out['mu_arr'], out['sigma_arr'], variables['uncertaintydata'], out['data'],shakemap)
    
    
    if plot_on == True:
        print('Plotting results')
        #plot(out, variables, voi, shakemap, stationdata)
        plot(out, variables, voi, shakemap,stationdata)
        
       
