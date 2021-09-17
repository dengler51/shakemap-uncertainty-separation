#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:32:38 2020

@author: dengler
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
sys.path.append('/Users/dengler/sm_corr_src')
sys.path.append('/Users/dengler/Documents/Codes/Loss with Uncertainty')
from setup import initialize
from loop import main
from realizations import realizations
from plotting import plot
from write_xml import write_grid_xml, readStation2

from pager_sc_realization import MC_pager_loss
from losspager.models.emploss import EmpiricalLoss
from losspager.models.exposure import Exposure
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
                                                                                                                                                                                                 
voi = 'PGA'
r = [15]
num_realizations = 0
cor_model = 'LB2013'
vscorr = True
plot_on = True



#np.random.seed(35626) 
np.random.seed(88532) 
#np.random.seed(332561) 
#np.random.seed(8824272) 

intra_flag = 1
station_flag = 0

#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/wenchuan_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/wenchuan_uncertainty.xml'
#station_file = 'none'

#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/events/nepal2_ps_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/events/nepal2_ps_uncertainty.xml'

#mmi_file = '/Users/dengler/shakemap_profiles/default/data/us20002926/Current/grid_alt.xml'
#unc_file = '/Users/dengler/shakemap_profiles/default/data/us20002926/Current/uncertainty_alt_phi.xml'
#tau_file = '/Users/dengler/shakemap_profiles/default/data/us20002926/Current/uncertainty_alt_tau.xml'
#

mmi_file = '/Users/dengler/shakemap_profiles/default/data/ci38457511/Current/grid_alt.xml'
unc_file = '/Users/dengler/shakemap_profiles/default/data/ci38457511/Current/uncertainty_alt_w.xml'
tau_file = '/Users/dengler/shakemap_profiles/default/data/ci38457511/Current/uncertainty_alt_b.xml'
#
#mmi_file = '/Users/dengler/shakemap_profiles/default/data/us2000ar20/current/grid_og.xml'
#unc_file = '/Users/dengler/shakemap_profiles/default/data/us2000ar20/current/uncertainty_og.xml'
#tau_file = '/Users/dengler/shakemap_profiles/default/data/us2000ar20/current/uncertainty_alt_tau.xml'

#mmi_file = '/Users/dengler/shakemap_profiles/default/data/hv70116556/current/Products/grid.xml'
#unc_file = '/Users/dengler/shakemap_profiles/default/data/hv70116556/current/Products/uncertainty.xml'
#
#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/us20002926_enhanced_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/us20002926_enhanced_uncertainty.xml'
#station_file = 'none'

#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/northridge_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/northridge_uncertainty.xml'
#station_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/northridge_stationlist.xml'

# mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/SMU_tests/EQ_data/JP/official20110311054624120_30/grid_nd.xml'
# unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/SMU_tests/EQ_data/JP/official20110311054624120_30/uncertainty.xml'
# tau_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/SMU_tests/EQ_data/JP/official20110311054624120_30/uncertainty_b_nd.xml'
#mmi_file = 'Inputs/grid.xml'
#unc_file = 'Inputs/uncertainty.xml'
#station_file = 'Inputs/stationlist2.xml'

#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/nr_nosta_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/nr_nosta_uncertainty.xml'
#station_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/nr_nosta_stationlist.xml'
#
#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/nr_sgmo_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/nr_sgmo_uncertainty.xml'
#station_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/nr_sgmo_stationlist.xml'

#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/lp_nosta_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/lp_nosta_uncertainty.xml'
#station_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/lp_nosta_stationlist.xml'
#
#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/lp_sgmo_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/lp_sgmo_uncertainty.xml'
#station_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/lp nr grids/lp_sgmo_stationlist.xml'
#
#mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/events/nepal2_ps_grid.xml'
#unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/events/nepal2_ps_uncertainty.xml'
#station_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/events/nr_sgmo_stationlist.xml'

# mmi_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/PCE Monte Carlo codes/pce_data/grid_og.xml'
# unc_file = '/Users/dengler/Documents/Codes/Loss with Uncertainty/PCE Monte Carlo codes/pce_data/uncertainty_og.xml'


#range_lat = [187,425]
#range_lon = [250,450]

#Gorkha
#range_lat = [60,415]
#range_lon = [100,460]

#Gorkha_alt
#range_lat = [67,422]
#range_lon = [107,467]

#Gorkha 2
#range_lat = [26,426]
#range_lon = [12,467]

#Gorkha alt 2
#range_lat = [32,432]
#range_lon = [19,474]

#range_lat = [50,310]
#range_lon = [50,365]
#Maule_zoom
#range_lat = [141,196]
#range_lon = [332,363]

#Maule
#range_lat = [84,408]
#range_lon = [60,462]

#range_lat = [100,525]
#range_lon = [100,567]
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
#        unc_INTRA = ShakeGrid(unc_file, variable= 'GMPE_INTRA_STD%s' % voi)
#        unc_INTER = ShakeGrid(unc_file, variable= 'GMPE_INTER_STD%s' % voi)
        unc_INTRA = ShakeGrid(unc_file, variable= 'STD%s' % voi)
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
#    f = open('rand.txt', 'r')
#    for i in range(0,variables['N']*variables['M']):
#        rand.append(float(f.readline()))
    
    rand = np.random.randn(variables['N']*variables['M'])
    
    # Truncate to +/- 3 sigma
    for i in range(0,variables['N']*variables['M']):
        while abs(rand[i])>3:
            rand[i]=np.random.randn(1)
#        if abs(rand[i]) > 2:
#            rand[i] = (4-np.sign(rand[i])*rand[i])
#    
    print('Calling main')
    out = main(variables, r, voi, rand, cor_model, vscorr)
    
    if num_realizations > 0:
        print('Computing realizations')
        realizations(num_realizations, radius, variables['N'], variables['M'], out['grid_arr'],
                 out['mu_arr'], out['sigma_arr'], variables['uncertaintydata'], out['data'],shakemap)
    
    #mmi_file = open('/Users/dengler/shakemap_profiles/default/data/ci38457511/Current/mmi_sta_locs.txt','r')
    #sta_file = open('/Users/dengler/shakemap_profiles/default/data/ci38457511/Current/inst_sta_locs.txt','r')

    if plot_on == True:
        print('Plotting results')
        #plot(out, variables, voi, shakemap, stationdata)
        plot(out, variables, voi, shakemap,stationdata)

mmi_file_name = 'grid_real.xml'        
        
write_grid_xml(shakemap,variables,voi,out['data_new'],mmi_file_name) 


fat_rate_file = '/Users/dengler/pager_src/losspager/data/fatality.xml'
mmi_MC_file = '/Users/dengler/sm_corr_src/'+mmi_file_name
pop_file = '/Users/dengler/pager_src/loss_pager/model_data/lspop2012.flt'
iso_file = '/Users/dengler/Documents/Codes/Loss With Uncertainty/isogrid.bil'

pop_year = int(str(shakemap.getAttributes()['event']['event_timestamp'])[0:4])        

empfat = EmpiricalLoss.fromXML(fat_rate_file)
expmodel = Exposure(pop_file,pop_year,iso_file,None)
    
        
# losses_MC,exposure_MC =MC_pager_loss(mmi_MC_file,expmodel,empfat) 
# losses_pager,exposure_pager =MC_pager_loss(mmi_file,expmodel,empfat) 

print(losses_MC)
print(' ')
print(exposure_MC)

#print(exposure)
#print(' ')
#print(' ')
#print(losses)

       
#xx = np.linspace(0,100,1000)
#yy = np.exp(-3*xx/25.7)
#plt.plot(xx,yy,'k',linewidth=3)
#plt.xlabel(r'Site Distance $\Delta$ (km)',fontsize = 16)
#plt.ylabel('Correlation', fontsize =16)
#plt.title('Correlation Function for MMI (PSA1.0) \n Using Jayaram & Baker (2009) Model',fontsize=18)
#plt.text(x=50,y=.8,s=r'$corr(\Delta)=e^{-\frac{3\Delta}{25.7}}$',fontsize=22)

#fig = plt.figure(figsize=(8,4))
#xx = np.linspace(0,100,1000)
#yy = .32*np.exp((-3.0/20)*xx)+.38*np.exp((-3.0/70)*xx) 
#plt.plot(xx,yy,'k',linewidth=3)
#plt.xlabel(r'Site Distance $\Delta$ (km)',fontsize = 16)
#plt.ylabel('Correlation', fontsize =16)
#plt.title('Correlation Function for MMI (PSA1.0) \n Using Loth & Baker (2013) Model',fontsize=18)
#plt.text(x=10,y=.6,s=r'$corr(\Delta,T=1s)=.32e^{-\frac{3\Delta}{20}}+.38e^{-\frac{3\Delta}{70}}+.30*I_{h=0}$',fontsize=16)