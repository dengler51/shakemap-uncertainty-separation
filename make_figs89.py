#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:09:06 2020

@author: dengler
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, uniform
import random
from numpy.linalg import cholesky
from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'inline')


# %%-------------------------------------------------------------------------
# Goda/Atkinson (2009) spatial correlation for PGA
# Input: a numpy array of distances in km
# Output: spatial correlation coefficients (from 0.0 to 1.0)
def ga_spatial_correlation(dists):
#    alpha = 0.060
#    beta = 0.283
#    gamma = 5.0
#    nal = -1.0 * alpha
#    gm1 = gamma - 1.0
#    cor = 1.0 - np.sqrt(1.0 - np.maximum(
#        gamma * np.exp(nal * np.power(dists, beta)) - gm1, 0))
    
    # This makes the function a simple exponential; a divisor of 
    # around 6 is similar to the Goda-Atkinson model (above),
    # larger values mean greater correlation, smaller values mean
    # shorter correlation
    cor = np.exp(-dists/6.0)
    #cor = np.exp(-dists/.10)
    return np.clip(cor, 0, 1)

def lb_spatial_correlation(dists,im1,im2):
    B1 = np.zeros_like(im1)
    B2 = np.zeros_like(im1)
    B3 = np.zeros_like(im1)
    for i in range(len(im1[:,0])):
        for j in range(len(im2[0,:])):
            if im1[i,j] == .01:
                if im2[i,j] == .01:
                    B1[i,j] =0.305
                    B2[i,j] =0.315
                    B3[i,j] =0.38
                elif im2[i,j] == 1:
                    B1[i,j] =0.16
                    B2[i,j] =0.17
                    B3[i,j] =0.04
            elif im1[i,j] == 1:
                if im2[i,j] == .01:
                    B1[i,j] =0.16
                    B2[i,j] =0.17
                    B3[i,j] =0.04
                elif im2[i,j] == 1:
                    B1[i,j] =0.32
                    B2[i,j] =0.38
                    B3[i,j] =0.30
                    
    afact = -3.0 / 20.0  # noqa
    bfact = -3.0 / 70.0  # noqa
    
     #Compute the correlation coefficient (Equation 42)
    
    corr = (B1)*np.exp(dists*afact) + (B2)*np.exp(dists*bfact) + (dists==0)*B3
    return np.clip(corr,0,1)


# %%-------------------------------------------------------------------------

# Notation:
# im1 = predictions
# im2 = observations
# phi = within-event standard deviation
# tau = between-event standard deviation
    
num_sims = 4000

# Number of points in the X_1 (prediction) vector
ngrid = 221
sgrid = 200

# Number of stations with data
nsta =14

#random.seed(15)

cutoff_dist = 120
cutoff_flag = 0

# Locations (in km) of the points in X_1
im1_dist = np.linspace(0.0, 220.0, ngrid).reshape((-1, 1))

# Select observation indices
#im2_dist = np.sort(np.array(random.sample(range(sgrid),nsta)))
#im2_dist = np.concatenate([ind_obs,np.array([300,330])])
#nsta +=2
#ind_obs = [2,450]
## Original toy problem 
im2_dist = np.array([10, 20, 30, 40, 44, 48, 52,
                 90, 98, 106, 114, 160, 168, 176], dtype=int)
nsta = len(im2_dist)


# Locations of stations

ind_obs = np.where(im1_dist == im2_dist)[0]

im2_dist = im1_dist[ind_obs]


# Select off-grid points (randomly selects a random number of stations and distributes
# them uniformly between the grid points)
num_off_grid_sta = random.sample(range(1,int(nsta/2)+1),1)[0]
#num_off_grid_sta = nsta
#off_grid_sta_ind = np.sort(np.array(random.sample(range(nsta),num_off_grid_sta)))
#im2_dist[off_grid_sta_ind] += uniform.rvs(size=(num_off_grid_sta,1))


off_grid_sta_ind = np.arange(nsta)[[-1,-2,-3]]
#off_grid_sta_ind = []
im2_dist[[off_grid_sta_ind]] += 0.5
#off_grid_sta_ind = []

#Indices of the stations which are on grid
on_grid_sta_ind = np.delete(np.arange(nsta),off_grid_sta_ind)



# "GMPE" estimates at the locations in the X_1 array; here we use
# a constant (or a linear ramp) for demonstration purposes
mu_im1 = np.full_like(im1_dist, -1.0)
# mu_im1 = -1.0 + 0.001 * im1_dist
#m1 = -(3.5)/(len(im1_dist)-100)
#mu_im1 = np.concatenate([7.5*np.ones(100),m1*np.arange(100,len(im1_dist)) +7.5-m1*100 ]).reshape([-1,1])

# "GMPE" phi and tau at locations in X_1; (increasing with distance for
# demonstration purposes)
#phi_im1 = 0.6 - 0.001*im1_dist
tau_im1 = 0.25 + 0.001*im1_dist
#tau_im1 = .25 -np.abs(0.001*(im1_dist-50))
phi_im1 = np.full_like(im1_dist, 0.6)
#tau_im1 = np.full_like(im1_dist, 0.15)
#tau_im1 = np.exp(-im1_dist/400)*(.6 + .1*np.sin(im1_dist/20 +10))
#phi_im1 = np.exp(-im1_dist/400)*(.6+.1*np.sin(im1_dist/40 +10))
#m1_tau = -(1.1-.7)/80
#m2_tau = .1/60
#m3_tau = .05/(len(im1_dist)-240)
#m1_phi = -(1.5-1.2)/50
#m3_phi = .05/(len(im1_dist)-240)
#tau_im1 = np.concatenate([1.1*np.ones(100), m1_tau*np.arange(100,180) +1.1-m1_tau*100,\
#                          m2_tau*np.arange(180,240) +.7-m2_tau*180, m3_tau*np.arange(240,len(im1_dist)) +.35-m3_tau*240]).reshape([-1,1])
#phi_im1 = np.concatenate([1.5*np.ones(100), m1_phi*np.arange(100,150) +1.5-m1_phi*100,\
#                          1.2*np.ones(90), m3_phi*np.arange(240,len(im1_dist)) +.8-m3_phi*240]).reshape([-1,1])

mu_im2 = mu_im1[ind_obs]
phi_im2 = phi_im1[ind_obs]
tau_im2 = tau_im1[ind_obs]

##Toy problem extra noise
obs_extra_sig = np.array([0.3, 0, 0, 0, 0, 0, 0,
                           0.5, 0.5, 0.5, 0.5, 0, 0, 0]).reshape(-1, 1)
noisy_sta_ind = np.arange(nsta)[obs_extra_sig[:,0]>0]

noisy_off_grid_inds = []
noisy_on_grid_inds = []
not_noisy_on_grid_inds = []
for i in range(nsta):
    if i in noisy_sta_ind:
        if i in off_grid_sta_ind:
            noisy_off_grid_inds.append(i)
        else:
            noisy_on_grid_inds.append(i)
    else:
        if i not in off_grid_sta_ind:
            not_noisy_on_grid_inds.append(i)
noisy_off_grid_inds = np.array(noisy_off_grid_inds)
noisy_on_grid_inds = np.array(noisy_on_grid_inds)   

not_noisy_on_grid_inds = np.array(ind_obs[not_noisy_on_grid_inds])         
        


#obs_extra_sig = uniform.rvs(size=(nsta,1))    
#imstr =['PGA','PGA','PGA','PGA','PGA','PGA','PGA',
#        'MMI','MMI','PGA','MMI','PGA','PGA','PGA'] 
imstr =['PGA','PGA','PGA','PGA','PGA','PGA','PGA',
        'PGA','PGA','PGA','PGA','PGA','PGA','PGA'] 
#imstr =['MMI','PGA','PGA','PGA','PGA','PGA','PGA',
#        'MMI','MMI','MMI','MMI','PGA','PGA','PGA']
im2_measure = np.zeros_like(obs_extra_sig)
for i in range(nsta):
    if imstr[i]=='PGA':
        im2_measure[i] = 0.01
    elif (imstr[i]=='MMI') | (imstr[i]=='PGV'):
        im2_measure[i] = 1

meas1_22 = im2_measure*np.ones_like(im2_measure.T)
meas2_22 = im2_measure.T*np.ones_like(im2_measure)
sigma22 = phi_im2 * phi_im2.T
dist22 = np.abs(im2_dist.T - im2_dist)
cov_dW2dW2 = lb_spatial_correlation(dist22,meas1_22,meas2_22)*sigma22
np.fill_diagonal(cov_dW2dW2, cov_dW2dW2.diagonal() + obs_extra_sig[:,0]**2)

#Draw the bias residiual and within residual realizations
#np.random.seed(2)
between_stand_norm_real = norm.rvs(size=1)
between_stand_norm_real = 1.5
#between_stand_norm_real = 2
between_res = norm.rvs(size=1)
between_res = 2
biased_mu = mu_im1 + tau_im1*between_stand_norm_real
between_res = between_stand_norm_real*tau_im2
#np.random.seed(11)
#within_res = cholesky(cov_dW2dW2).dot(norm.rvs(size=(nsta,1)))



# Observations
#im2 = mu_im2 + between_res + within_res 

#Original toy problem data
im2 =-.1+  np.array([-.8, -1.2, -0.9, -0.9, -1.1, -0.8, -1.2,
                 -.7, -0.9, -0.8, -0.9, -1.1, -1.2, -1.3]).reshape((-1, 1))


   
# Additional sigma for some of the observations
#obs_extra_sig = np.full_like(phi_im2,0.0) #No noise

#num_obs_extra_sig = random.sample(range(1,int(nsta/2)+1),1)[0]
#noisy_sta_ind = np.sort(np.array(random.sample(range(nsta),num_obs_extra_sig)))
#obs_extra_sig[noisy_sta_ind] = 0.5




# The raw residuals
zeta = im2 - mu_im2


#%%
# %%-------------------------------------------------------------------------
# Do the bias correction


# Now add the extra uncertainty
cov_dW2dW2_inv = np.linalg.inv(cov_dW2dW2)


# The event term and its variance (for the observation points)
if cutoff_flag == 1:
    dix = (im2_dist <= cutoff_dist).ravel()
else:
    dix = np.ones(len(im2_dist))==1
    
elems = np.delete(np.arange(len(im2_dist)),np.arange(len(im2_dist))[dix])
cov_inv_dix = np.delete(np.delete(cov_dW2dW2_inv,elems,axis=0),elems,axis=1)  
tau_im2_dl = tau_im2[dix].reshape([-1,1])   
zeta_dl = zeta[dix].reshape([-1,1])

bias_denom = (tau_im2_dl).T.dot(cov_inv_dix.dot(tau_im2_dl))   
bias_num = (tau_im2_dl).T.dot(cov_inv_dix.dot(zeta_dl))

sig_delta_B_im2 =tau_im2*np.sqrt(1.0/(1.0 + bias_denom) )
delta_B_im2 = 1.0/tau_im2*sig_delta_B_im2**2*bias_num


# The event term, its variance and the total within-event standard
# deviation (for the predictions)
sig_delta_B_im1 =tau_im1*np.sqrt(1.0/(1.0 + bias_denom) )
delta_B_im1 = 1.0/tau_im1*sig_delta_B_im1**2*bias_num

# %%-------------------------------------------------------------------------
# Solve for the conditional mean and covariance

# Distance matrices
dist11 = np.abs(im1_dist.T - im1_dist)
dist12 = np.abs(im1_dist - im2_dist.T)

cov_dW1dW1 = lb_spatial_correlation(dist11,0.01*np.ones_like(dist11),0.01*np.ones_like(dist11)) * (phi_im1 * phi_im1.T)
cov_dW1dW2 = lb_spatial_correlation(dist12,0.01*np.ones_like(dist12),im2_measure.T*np.ones_like(im1_dist)) * (phi_im1 * phi_im2.T) 

cov_dB1dB1 = tau_im1 * tau_im1.T
cov_dB1dB2 = tau_im1 * tau_im2.T
cov_dB2dB2 = tau_im2 * tau_im2.T

cov_dB1dB1_im2 = sig_delta_B_im1 * sig_delta_B_im1.T
cov_dB1dB2_im2 = sig_delta_B_im1 * sig_delta_B_im2.T
cov_dB2dB2_im2 = sig_delta_B_im2 * sig_delta_B_im2.T

# Covariance for im1 x im1 and im1 x im2
cov_im1im1 = cov_dW1dW1 + cov_dB1dB1
cov_im1im2 = cov_dW1dW2 + cov_dB1dB2 

# Build the updated (by bias) covariance matrix of the residuals and its inverse
cov_im2im2 = cov_dW2dW2 + cov_dB2dB2


# Now add the extra uncertainty
np.fill_diagonal(cov_im2im2, cov_im2im2.diagonal() + obs_extra_sig[:,0]**2)

# Invert
cov_im2im2_inv = np.linalg.inv(cov_im2im2)

# Regression coefficient matrix
rcmatrix_im = cov_im1im2.dot(cov_im2im2_inv)
rcmatrix_dW = cov_dW1dW2.dot(cov_dW2dW2_inv)


# Conditional normalized mean
mu_im1_im2 = mu_im1 + rcmatrix_im.dot(zeta)
mu_im1_im2 = mu_im1 + delta_B_im1 + rcmatrix_dW.dot(zeta - delta_B_im2)

# Conditional normalized covariance
#cov_im1im1_im2 = cov_im1im1 - rcmatrix_im.dot(cov_im1im2.T)
cov_im1im1_im2 = cov_dW1dW1 - rcmatrix_dW.dot(cov_dW1dW2.T) +\
               rcmatrix_dW.dot(cov_dB2dB2_im2.dot(rcmatrix_dW.T)) - \
               rcmatrix_dW.dot(cov_dB1dB2_im2.T) - cov_dB1dB2_im2.dot(rcmatrix_dW.T) +\
               cov_dB1dB1_im2
cov_im1im1_im2_within = cov_dW1dW1 - rcmatrix_dW.dot(cov_dW1dW2.T) #+ \
      #          cov_dB1dB1_im2
cov_im1im1_im2_between = rcmatrix_dW.dot(cov_dB2dB2_im2.dot(rcmatrix_dW.T)) - \
               rcmatrix_dW.dot(cov_dB1dB2_im2.T) - cov_dB1dB2_im2.dot(rcmatrix_dW.T) +\
               cov_dB1dB1_im2


# Conditional standard deviation (sqrt of the diagonal of the
# conditional covariance matrix)
sig_im1im1_im2_diag = np.sqrt(np.clip(np.diag(cov_im1im1_im2), 0, np.inf))
sig_im1im1_im2_diag_within = np.sqrt(np.clip(np.diag(cov_im1im1_im2_within), 0, np.inf))
sig_im1im1_im2_diag_between = np.sqrt(np.clip(np.diag(cov_im1im1_im2_between),0,np.inf))




#%% Monte Carlo Simulation
num_runs = 1
z_stds = np.zeros(num_runs)
for run in range(num_runs):


    gm_realizations = np.zeros((num_sims,len(im1_dist)))
    updated_cov_mat = np.delete(np.delete(cov_im1im1_im2_within,not_noisy_on_grid_inds,axis=0),not_noisy_on_grid_inds,axis=1)
    gm_realizations[:,not_noisy_on_grid_inds] = np.tile(mu_im1_im2[not_noisy_on_grid_inds].T,(num_sims,1))
    
    sig_im1im1_im2_diag_within[sig_im1im1_im2_diag<1E-8] = 0
    gm_realizations_alt = np.zeros((num_sims,len(im1_dist)))
    cov_im1im1_im2_within_alt = lb_spatial_correlation(dist11,0.01*np.ones_like(dist11),0.01*np.ones_like(dist11)) * (sig_im1im1_im2_diag_within[:,np.newaxis] * sig_im1im1_im2_diag_within[np.newaxis,:])
    updated_cov_mat_alt = np.delete(np.delete(cov_im1im1_im2_within_alt,not_noisy_on_grid_inds,axis=0),not_noisy_on_grid_inds,axis=1)
    gm_realizations_alt[:,not_noisy_on_grid_inds] = np.tile(mu_im1_im2[not_noisy_on_grid_inds].T,(num_sims,1))
    
    
    #Simulation random fields
    within_rands = norm.rvs(size=(len(im1_dist)-len(not_noisy_on_grid_inds),num_sims))
    
    between_event_residuals = sig_im1im1_im2_diag_between[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)][np.newaxis,:]*norm.rvs(size=(num_sims,1))
    within_event_residuals = np.dot(cholesky(updated_cov_mat) ,within_rands).T
    within_event_residuals_alt = np.dot(cholesky(updated_cov_mat_alt) ,within_rands).T
    
    
    
    mean_between = np.mean(between_event_residuals,axis=0)
    std_between = np.std(between_event_residuals,axis=0)
    
    mean_between_p = np.zeros(len(im1_dist))
    mean_between_p[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)]=mean_between
    std_between_p = np.zeros(len(im1_dist))
    std_between_p[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)]=std_between
    
    #between_res = np.zeros(len(im1_dist))
    within_res = np.zeros(len(im1_dist))
    within_res_alt = np.zeros(len(im1_dist))
    #within_std = np.zeros(len(im1_dist))
    #within_std[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)] = 
    #between_res[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)] = between_event_residuals[0,:]
    within_res[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)] = within_event_residuals[0,:]
    within_res_alt[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)] = within_event_residuals_alt[0,:]
    
    between_res = sig_im1im1_im2_diag_between*.5
    
    
    mean_within = np.mean(within_event_residuals,axis=0)
    std_within = np.std(within_event_residuals,axis=0)
    
    mean_within_p = np.zeros(len(im1_dist))
    mean_within_p[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)]=mean_within
    std_within_p = np.zeros(len(im1_dist))
    std_within_p[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)]=std_within
    
    
    gm_realizations[:,np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)] = mu_im1_im2[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)].T  + within_event_residuals+ between_event_residuals
    gm_realizations_alt[:,np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)] = mu_im1_im2[np.delete(np.arange(len(im1_dist)),not_noisy_on_grid_inds)].T  + within_event_residuals_alt+ between_event_residuals
    
    
    
    mean_gm = np.mean(gm_realizations,axis=0)
    std_gm = np.std(gm_realizations,axis=0)
    cov_gm = np.cov(gm_realizations.T)
    
    mean_gm_alt = np.mean(gm_realizations_alt,axis=0)
    std_gm_alt= np.std(gm_realizations_alt,axis=0)
    cov_gm_alt= np.cov(gm_realizations_alt.T)
    
    
    
    con_mean = np.mean(gm_realizations[:,300:],axis=1)
    con_std = np.std(gm_realizations[:,300:],axis=1)
    
    gm_no_data_realizations = mu_im1.T + np.dot(cholesky(ga_spatial_correlation(dist11)*(phi_im1 * phi_im1.T)+tau_im1 * tau_im1.T),norm.rvs(size=(len(phi_im1),num_sims))).T
    #gm_no_data_realizations = mu_im1.T + np.dot(cholesky(ga_spatial_correlation(dist11)*(phi_im1 * phi_im1.T)),norm.rvs(size=(len(phi_im1),num_sims))).T
    
    
    
    mean_gm_no_data = np.mean(gm_no_data_realizations,axis=0)
    std_gm_no_data = np.std(gm_no_data_realizations,axis=0)
    
    con_mean_no_data = np.mean(gm_no_data_realizations[:,300:],axis=1)
    con_std_no_data = np.std(gm_no_data_realizations[:,300:],axis=1)
    
    # #%%
    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot()
    ax1.plot(im1_dist, mu_im1,'k:', label=r'GMPE $\mu_{Y_1}$',linewidth=2.5)
    ax1.plot(im1_dist, mu_im1_im2,'b-',alpha=.8, label=r'$\mu_{Y_1|y_2}$',linewidth=2.5,zorder=20)
    ax1.plot(im1_dist, gm_realizations[0,:],'k-',alpha=.3, label='Random field')
    ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+sig_im1im1_im2_diag,mu_im1_im2[:,0]-sig_im1im1_im2_diag,facecolor='black',alpha=.2,label='+/- 1 STD region')
    ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+2*sig_im1im1_im2_diag,mu_im1_im2[:,0]+sig_im1im1_im2_diag,facecolor='black',alpha=.1,label='+/- 2 STD region')
    ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]-sig_im1im1_im2_diag,mu_im1_im2[:,0]-2*sig_im1im1_im2_diag,facecolor='black',alpha=.1)
    #ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
    #ax1.plot(im1_dist,mean_gm,label='Mean GM')
    
    #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
    #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')
    
    #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
    ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
    ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
    if noisy_on_grid_inds != []:
        ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
    if noisy_off_grid_inds != []:
        ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)
    
    
    #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
    #ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
    ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))    
     
    ax1.set_ylabel('log PGA', fontsize=16)
    ax1.set_xlabel('Location (km)',fontsize=16)
    ax1.tick_params(labelsize=16)
    #plt.ylim([-2,.5])
    
    
    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot()
    ax1.plot(im1_dist, mu_im1,'k:', label=r'GMPE $\mu_{Y_1}$',linewidth=2.5)
    ax1.plot(im1_dist, mu_im1_im2,'b-',alpha=.8, label=r'$\mu_{Y_1|y_2}$',linewidth=2.5,zorder=20)
    ax1.plot(im1_dist, gm_realizations[0,:].T,'k-',alpha=.35, label='Random fields')
    ax1.plot(im1_dist, gm_realizations[1:20,:].T,'k-',alpha=.35)
    ax1.plot(im1_dist,np.mean(gm_realizations[0:20,:],axis=0),'r-',label= 'Sample Mean',linewidth=2.5)
    ax1.plot(im1_dist, mu_im1_im2.squeeze() + np.std(gm_realizations[0:20,:],axis=0),'g-',label='+/- 1 Sample STD',linewidth=2.5)
    ax1.plot(im1_dist, mu_im1_im2.squeeze() - np.std(gm_realizations[0:20,:],axis=0),'g-',linewidth=2.5)

    ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+sig_im1im1_im2_diag,mu_im1_im2[:,0]-sig_im1im1_im2_diag,facecolor='black',alpha=.2,label='+/- 1 STD region')
    #ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
    #ax1.plot(im1_dist,mean_gm,label='Mean GM')
    
    #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
    #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')
    
    #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
    ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
    ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
    if noisy_on_grid_inds != []:
        ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
    if noisy_off_grid_inds != []:
        ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)
    
    
    #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
    #ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
    ax1.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.001, 1))  
     
    ax1.set_ylabel('log PGA', fontsize=16)
    ax1.set_xlabel('Location (km)',fontsize=16)
    ax1.tick_params(labelsize=16)
    plt.ylim([-2.5,.7])
    plt.xlim([5,35])
    
    # fig_zoom = plt.figure(figsize=(10, 6))
    # ax1 = fig_zoom.add_subplot()
    # ax1.plot(im1_dist, mu_im1,'k:', label=r'GMPE $\mu_{Y_1}$',linewidth=2.5)
    # ax1.plot(im1_dist, mu_im1_im2,'b-',alpha=.8, label=r'$\mu_{Y_1|y_2}$',linewidth=3,zorder=20)
    # ax1.plot(im1_dist, gm_realizations[0,:],'k-',alpha=.3, label='Random fields')
    # ax1.plot(im1_dist, gm_realizations[1:np.min([num_sims,20]),:].T,'k-',alpha=.3)
    # ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+sig_im1im1_im2_diag,mu_im1_im2[:,0]-sig_im1im1_im2_diag,facecolor='black',alpha=.1,label='+/- 1 STD region')
    
    # #ax1.plot(im1_dist,mean_gm,label='Mean GM')
    
    # #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
    # #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')
    
    # #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
    # ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 17,label=r'On Grid $y_2$',zorder=30)
    # #ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
    # if noisy_on_grid_inds != []:
    #     ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=15,label=r'Noisy $y_2$',zorder=31)
    # if noisy_off_grid_inds != []:
    #     ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)
    
    
    # #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
    # ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))    
    # ax1.set_ylabel('ln(PGA)', fontsize=16)
    # ax1.set_xlabel('Location (km)',fontsize=16)
    # ax1.tick_params(labelsize=16)
    # plt.xlim([7,33])
    # plt.ylim([-2.5,.5])
    
    
    
    
    
    fig = plt.figure(figsize=(10, 12))
    ax1 = fig.add_subplot(211)
    ax1.plot(im1_dist, mu_im1,'k:',linewidth=2.5, label=r'GMPE $\mu_{Y_1}$')
    ax1.plot(im1_dist, mu_im1_im2,'b--',alpha=.7, label=r'$\mu_{Y_1|y_2}$',linewidth=2.5,zorder=20)
    #ax1.plot(im1_dist, gm_realizations[0,:],':', label='Random fields')
    #ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
    ax1.plot(im1_dist,mean_gm,'r-',alpha=.4,linewidth=2.5,label=r'Sample $\mu_{Y_1|y_2}$')
    #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
    #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')
    
    #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
    ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
    ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
    if noisy_on_grid_inds != []:
        ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
    if noisy_off_grid_inds != []:
        ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)
    
    
    #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
    ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))    
    ax1.set_ylabel('ln(PGA)', fontsize=16)
    
    ax1.tick_params(labelsize=16)
    #ax1.set_ylim((-1.32, -0.25))
    # ax1.set_ylim((-2.32, -0.65))
    
    # Plot the standard deviation
    
    ax2 = fig.add_subplot(212)
    #ax2.plot(im1_dist, np.sqrt(phi_im1**2 + tau_im1**2), label='Prior total')
    
    #ax2.plot(im1_dist, sig_im1im1_im2_diag_og, label='Cond Total Std (Current Model)')
    ax2.plot(im1_dist, sig_im1im1_im2_diag,'b--',alpha=.7, label=r'$\sigma_{Y_1|y_2}$',linewidth=2.5,zorder=20)
    ax2.plot(im1_dist,std_gm,'r-',linewidth=2.5,alpha=.4,label=r'Sample $\sigma_{Y_1|y_2}$')
    #ax2.plot(im1_dist, sig_im1im1_im2_diag_within, label='Conditional Within Std',linewidth=2,zorder=21)
    #ax2.plot(im1_dist,sig_delta_B_im1,label='Updated bias STD')
    #ax2.plot(im1_dist, sig_im1im1_im2_diag_between, label='Conditional Between Std',linewidth=2,zorder=22)
    ax2.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$',zorder=32)
    ax2.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=8,zorder=33)
    if noisy_on_grid_inds != []:
        ax2.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $\sigma_{y_2}$',zorder=30)
    if noisy_off_grid_inds != []:
        ax2.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',zorder=31)
    
    
    #ax2.plot(im1_dist,std_within_p,label='Within STD',zorder = 28)
    #ax2.plot(im1_dist,std_between_p,label='Between STD',zorder=27)
    
    #ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
    ax2.legend(loc='upper right', bbox_to_anchor=(1.232, 1))
    
    ax2.set_xlabel('Location (km)', fontsize=16)
    ax2.set_ylabel('Standard deviation', fontsize=16)
    ax2.tick_params(labelsize=16)
    #ax2.set_ylim((-0.01, 0.75))
    
    # fig3 = plt.figure(figsize=(8,6))
    # plt.plot(im1_dist,mean_between_p,label='Between Mean',zorder=20)
    # plt.plot(im1_dist,mean_within_p,label='Within Mean')
    # plt.plot(im1_dist,mean_gm-mu_im1_im2[:,0],label='Sample Mean-Actual Mean')
    # #plt.plot(im1_dist,std_gm-sig_im1im1_im2_diag,label='Sample STD-Actual STD')
    # plt.legend()
    
    fig.savefig('/Users/dengler/Documents/uncertainty paper figures/Figure10.jpg', format='jpg', dpi=1200, bbox_inches="tight")

    mean_std_300 = np.mean(std_gm[300:])
    std_std_300 = np.std(std_gm[300:])
    
    mean_mean_300 = np.mean(mean_gm[300:])
    std_mean_300 = np.std(mean_gm[300:])
    
    mean_std_actual = np.mean(sig_im1im1_im2_diag[300:])
    std_std_actual = np.std(sig_im1im1_im2_diag[300:])
    
    mean_mean_actual = np.mean(mu_im1_im2[300:])
    std_mean_actual = np.std(mu_im1_im2[300:])
    
    con_mean_mean = np.mean(con_mean)
    con_mean_std = np.std(con_mean)
    
    con_std_mean = np.mean(con_std)
    con_std_std = np.std(con_std)
    
    betw_mean_std = np.mean(sig_im1im1_im2_diag_between[300:])
    
    mean_std_no_data_300 = np.mean(std_gm_no_data[300:])
    std_std_no_data_300 = np.std(std_gm_no_data[300:])
    
    mean_mean_no_data_300 = np.mean(mean_gm_no_data[300:])
    std_mean_no_data_300 = np.std(mean_gm_no_data[300:])
    
    
    con_mean_mean_no_data = np.mean(con_mean_no_data)
    con_mean_std_no_data = np.std(con_mean_no_data)
    
    con_std_mean_no_data = np.mean(con_std_no_data)
    con_std_std_no_data = np.std(con_std_no_data)
    
    
    z = (mean_gm - mu_im1_im2.squeeze())/sig_im1im1_im2_diag.squeeze()
    z = np.delete(z,np.isnan(z))
    z = np.delete(z,np.isinf(z))
    z_stds[run] = np.std(z)
   # print(np.std(z))
    
  #  zs = np.array([.0134,.0145,0.01268 ])
# fig = plt.figure(figsize=(10, 6))
# ax1 = fig.add_subplot()
# ax1.plot(im1_dist, mu_im1,'k:', label=r'GMPE $\mu_{Y_1}$',linewidth=2.5)
# ax1.plot(im1_dist, mu_im1_im2,'b-',alpha=.8, label=r'$\mu_{Y_1|y_2}$',linewidth=2.5,zorder=20)
# ax1.plot(im1_dist, gm_realizations[0,:],'k-',alpha=.3, label='Random field')
# ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+sig_im1im1_im2_diag,mu_im1_im2[:,0]-sig_im1im1_im2_diag,facecolor='black',alpha=.1,label='+/- 1 STD region')
# #ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
# #ax1.plot(im1_dist,mean_gm,label='Mean GM')

# #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
# #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')

# #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
# ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
# ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
# if noisy_on_grid_inds != []:
#     ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
# if noisy_off_grid_inds != []:
#     ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)


# #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
# #ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
# ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))    
 
# ax1.set_ylabel('ln(PGA)', fontsize=16)
# ax1.set_xlabel('Location (km)',fontsize=16)
# ax1.tick_params(labelsize=16)
# plt.ylim([-2,.5])


# #fig = plt.figure(figsize=(8, 5))
# #ax1 = fig.add_subplot()
# #ax1.plot(im1_dist, mu_im1,'k:', label=r'GMPE $\mu_{Y_1}$',linewidth=2.5)
# #ax1.plot(im1_dist, mu_im1_im2,'k-',alpha=.7, label=r'$\mu_{Y_1|y_2}$',linewidth=2.5,zorder=20)
# #ax1.plot(im1_dist, mu_im1[:,0] + between_res + within_res,'k-',alpha=.3, label='Random field')
# #ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+sig_im1im1_im2_diag,mu_im1_im2[:,0]-sig_im1im1_im2_diag,facecolor='black',alpha=.1,label='+/- 1 Total\nSTD region')
# ##ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
# ##ax1.plot(im1_dist,mean_gm,label='Mean GM')
# #
# ##ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
# ##ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')
# #
# ##ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
# #ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
# #ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
# #if noisy_on_grid_inds != []:
# #    ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
# #if noisy_off_grid_inds != []:
# #    ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)
# #
# #
# ##ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
# ##ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
# #ax1.legend(loc='upper right',fontsize = 14, bbox_to_anchor=(1.35, 1))    
# # 
# #ax1.set_ylabel('log PGA', fontsize=16)
# #ax1.set_xlabel('Location (km)',fontsize=16)
# #ax1.tick_params(labelsize=16)
# #plt.ylim([-2,.5])
# #plt.xlim([-1,1+im1_dist[-1,0]])
# #
# #fig = plt.figure(figsize=(8, 5))
# #ax1 = fig.add_subplot()
# #ax1.plot(im1_dist,np.zeros(len(im1_dist)),'k:',linewidth=2.6)
# #ax1.plot(im1_dist, within_res,'r-',alpha=.6, label='Within \nRealization',linewidth=2.8)
# #ax1.fill_between(im1_dist[:,0],sig_im1im1_im2_diag_within,-sig_im1im1_im2_diag_within,facecolor='red',alpha=.3,label='+/- 1 Within\nSTD region')
# #ax1.plot(im1_dist, between_res,'b-',alpha=.6, label='Between \nRealization',linewidth=2.8)
# #ax1.fill_between(im1_dist[:,0],sig_im1im1_im2_diag_between,-sig_im1im1_im2_diag_between,facecolor='blue',alpha=.3,label='+/- 1 Between\nSTD region')
# ##ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
# ##ax1.plot(im1_dist,mean_gm,label='Mean GM')
# #
# ##ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
# ##ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')
# #
# ##ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
# ##ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
# ##ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
# ##if noisy_on_grid_inds != []:
# ##    ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
# ##if noisy_off_grid_inds != []:
# ##    ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)
# #
# #
# ##ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
# ##ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
# #ax1.legend(loc='upper right',fontsize = 14, bbox_to_anchor=(1.366, 1))    
# # 
# #ax1.set_ylabel('log PGA', fontsize=16)
# ##ax1.set_xlabel('Location (km)',fontsize=16)
# #ax1.tick_params(labelsize=16)
# #plt.ylim([-1.1,1.3])
# #plt.xlim([-1,1+im1_dist[-1,0]])

# ###########################################
# fig = plt.figure(figsize=(18, 6.5))
# ax1 = fig.add_subplot()
# ax1.plot(im1_dist, mu_im1,'k:', label=r'GMPE $\mu_{Y_1}$',linewidth=2.5)
# ax1.plot(im1_dist, mu_im1_im2,'k-',alpha=1, label=r'$\mu_{Y_1|y_2}$',linewidth=2.5,zorder=20)
# ax1.plot(im1_dist, mu_im1_im2[:,0]+ between_res ,'b-',alpha=.7, label='Mean + Between',linewidth=2.8)
# ax1.plot(im1_dist, mu_im1_im2[:,0] + between_res + within_res,'m-',alpha=.7,linewidth=2.8, label='Random field')
# #ax1.fill_between(im1_dist[:,0],mu_im1_im2[:,0]+sig_im1im1_im2_diag,mu_im1_im2[:,0]-sig_im1im1_im2_diag,facecolor='black',alpha=.1,label='+/- 1 Total\nSTD region')
# #ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
# #ax1.plot(im1_dist,mean_gm,label='Mean GM')

# #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
# #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')

# #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
# ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
# ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
# if noisy_on_grid_inds != []:
#     ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
# if noisy_off_grid_inds != []:
#     ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)


# #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
# #ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
# ax1.legend(loc='upper left',fontsize = 18, bbox_to_anchor=(1., 1),frameon=False)    
 
# ax1.set_ylabel('log PGA', fontsize=24)
# ax1.set_xlabel('Location (km)',fontsize=24)
# ax1.tick_params(labelsize=20)
# plt.ylim([-2,.5])
# plt.xlim([-1,1+im1_dist[-1,0]])

# fig = plt.figure(figsize=(18, 6.5))
# ax1 = fig.add_subplot()
# ax1.plot(im1_dist,np.zeros(len(im1_dist)),'k:',linewidth=2.6)
# ax1.fill_between(im1_dist[:,0],sig_im1im1_im2_diag_between,-sig_im1im1_im2_diag_between,facecolor='blue',alpha=.3,label='+/- 1 Between\nSTD region',zorder = 3)
# ax1.fill_between(im1_dist[:,0],sig_im1im1_im2_diag_within,-sig_im1im1_im2_diag_within,facecolor='red',alpha=.3,label='+/- 1 Within\nSTD region',zorder = 2)
# ax1.plot(im1_dist, between_res,'b-',alpha=1, label='Between \nRealization',linewidth=2.8,zorder=4)
# ax1.plot(im1_dist, within_res,'r-',alpha=.6, label='Within \nRealization',linewidth=2.8)

# #ax1.plot(im1_dist, gm_realizations[1:,:].T,':')
# #ax1.plot(im1_dist,mean_gm,label='Mean GM')

# #ax1.plot(im1_dist, mu_im1+delta_B_im1, label='Updated biased mu')
# #ax1.plot(im1_dist, mu_im1+delta_B_im1_og, label='Updated bias (Current)')

# #ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
# #ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$',zorder=30)
# #ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk', label=r'Off Grid $y_2$',markersize=8,zorder=32)
# #if noisy_on_grid_inds != []:
# #    ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$',zorder=31)
# #if noisy_off_grid_inds != []:
# #    ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red',zorder=33)


# #ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
# #ax1.legend(loc='upper right', bbox_to_anchor=(1.232, 1))   
# ax1.legend(loc='upper left',fontsize = 18, bbox_to_anchor=(1, 1),frameon=False)    
 
# ax1.set_ylabel('log PGA', fontsize=24)
# ax1.set_xlabel('Location (km)',fontsize=24)
# #ax1.set_xlabel('Location (km)',fontsize=16)
# ax1.tick_params(labelsize=20)
# plt.ylim([-1.2,1.3])
# plt.xlim([-1,1+im1_dist[-1,0]])
