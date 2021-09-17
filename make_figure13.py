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
    nug = .0
    cor = (1-nug)*np.exp(-dists/6.0) + nug*(dists==0)
    return np.clip(cor, 0, 1)

def lb_spatial_correlation(dists,im1,im2):
    B1 = np.zeros_like(im1)
    B2 = np.zeros_like(im1)
    B3 = np.zeros_like(im1)
    for i in range(len(im1[:,0])):
        for j in range(len(im2[0,:])):
            if im1[i,j] == .01:
                if im2[i,j] == .01:
                    B1[i,j] =0.29
                    B2[i,j] =0.47
                    B3[i,j] =0.24
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

def gmm_model(M,R):
    return -0.152+0.859*M - 1.803*np.log(R+25)


         
# %%-------------------------------------------------------------------------

# Notation:
# im1 = predictions
# im2 = observations
# phi = within-event standard deviation
# tau = between-event standard deviation

mag = 7.8

# Number of points in the X_1 (prediction) vector
ngrid = 261
sgrid = 200

# Number of stations with data
nsta = 14
#nsta = 22

#random.seed(15)

cutoff_dist = 100
cutoff_flag =0

# Locations (in km) of the points in X_1
im1_dist = np.linspace(0.0, 261.0, ngrid).reshape((-1, 1))

# Select observation indices
ind_obs = np.sort(np.array(random.sample(range(sgrid),nsta)))
#ind_obs = np.concatenate([ind_obs,np.array([300,330])])
#nsta +=2
#ind_obs = [2,450]
## Original toy problem 
ind_obs = np.array([10, 20, 30, 40, 44, 48, 52,
                  90, 98, 106, 114, 160, 168, 176], dtype=int)
    
#ind_obs = np.array([10, 20, 30, 40, 44, 48, 52,
#                 90,90,90,90,90,90,90,90,90, 98, 106, 114, 160, 168, 176], dtype=int)    
#ind_obs = np.array([10, 20, 30, 40, 44, 48, 52,
#                 90,92,93,94,96, 98, 102,104,106, 107,111,114, 160, 168, 176], dtype=int)
    

nsta = len(ind_obs)

# Locations of stations
im2_dist = im1_dist[ind_obs]


# Select off-grid points (randomly selects a random number of stations and distributes
# them uniformly between the grid points)
num_off_grid_sta = random.sample(range(1,int(nsta/2)+1),1)[0]
#num_off_grid_sta = nsta
# off_grid_sta_ind = np.sort(np.array(random.sample(range(nsta),num_off_grid_sta)))
# im2_dist[off_grid_sta_ind] += uniform.rvs(size=(num_off_grid_sta,1))

off_grid_sta_ind = np.arange(nsta)[-3:]
im2_dist[[off_grid_sta_ind]] += 0.5
#im2_dist[:] +=.5
#off_grid_sta_ind = np.arange(nsta)
##im2_dist[[1,2,3,4,5,6,11,12,13]] +=.5
#im2_dist[[0,7,8,9,10]]+=.1
#off_grid_sta_ind = np.array([1,2,3,4,5,6,11,12,13])

#Indices of the stations which are on grid
on_grid_sta_ind = np.delete(np.arange(nsta),off_grid_sta_ind)



# "GMPE" estimates at the locations in the X_1 array; here we use
# a constant (or a linear ramp) for demonstration purposes
mu_im1 = np.full_like(im1_dist, -1.0)
#mu_im1 = gmm_model(mag,im1_dist)

#mu_im1 = -1 - 0.001 * im1_dist
#m1 = -(3.5)/(len(im1_dist)-100)
#mu_im1 = np.concatenate([7.5*np.ones(100),m1*np.arange(100,len(im1_dist)) +7.5-m1*100 ]).reshape([-1,1])

# "GMPE" phi and tau at locations in X_1; (increasing with distance for
# demonstration purposes)
phi_im1 = 0.4 + 0.001*im1_dist
tau_im1 = 0.25 + 0.001*im1_dist
#tau_im1 = .6 -np.abs(0.001*(im1_dist-50))
phi_im1 = np.full_like(im1_dist, 0.6)
tau_im1 = np.full_like(im1_dist, 0.3)
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

# Additional sigma for some of the observations
obs_extra_sig = np.full_like(phi_im2,0.0) #No noise
#
# num_obs_extra_sig = random.sample(range(1,int(nsta/2)+1),1)[0]
# noisy_sta_ind = np.sort(np.array(random.sample(range(nsta),num_obs_extra_sig)))
# obs_extra_sig[noisy_sta_ind] = .5
#
###Toy problem extra noise
obs_extra_sig = np.array([0.3, 0, 0, 0, 0, 0, 0,
                            0.5, 0.5, 0.5, 0.5, 0, 0, 0]).reshape(-1, 1)
##obs_extra_sig = np.array([0.3, 0, 0, 0, 0, 0, 0,
#                           0.2,.5,.2,.5,.2,.2,.5,.2,.5, 0.5, 0.5, 0.5, 0, 0, 0]).reshape(-1, 1)
    
#obs_extra_sig = np.array([0.3, 0, 0, 0, 0, 0, 0,
#                           0.2, 0.2, 0.2, 0.3, 0.2, 0.4, 0.5, 0.2, 0.2, 0.3, 0.4, 0.5, 0, 0, 0]).reshape(-1, 1)

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
#im2_measure = .01*np.ones_like(obs_extra_sig)
#im2_measure[(imstr=='MMI') | (imstr == 'PGV')] = 1
#im2_measure = .01*np.ones_like(obs_extra_sig)
#im2_measure[obs_extra_sig==.5] = 1


meas1_22 = im2_measure*np.ones_like(im2_measure.T)
meas2_22 = im2_measure.T*np.ones_like(im2_measure)
sigma22 = phi_im2 * phi_im2.T
dist22 = np.abs(im2_dist.T - im2_dist)
cov_dW2dW2 = lb_spatial_correlation(dist22,meas1_22,meas2_22)*sigma22
np.fill_diagonal(cov_dW2dW2, cov_dW2dW2.diagonal() + obs_extra_sig[:,0]**2)

#Draw the bias residiual and within residual realizations
#np.random.seed(2)
between_stand_norm_real = norm.rvs(size=1)
#between_stand_norm_real = 2
between_res = norm.rvs(size=1)
biased_mu = mu_im1 + tau_im1*between_stand_norm_real
between_res = between_stand_norm_real*tau_im2
#np.random.seed(11)
within_res = cholesky(cov_dW2dW2).dot(norm.rvs(size=(nsta,1)))



# Observations
# im2 = mu_im2 + between_res + within_res 

#Original toy problem data
im2 =-.1+  np.array([-.8, -1.2, -0.9, -0.9, -1.1, -0.8, -1.2,
                  -.7, -0.9, -0.8, -0.9, -1.1, -1.2, -1.3]).reshape((-1, 1))
    
#im2 =  np.array([-0.8, -1.2, -0.9, -0.9, -1.1, -0.8, -1.2,
#                 -0.7,-.8,-.7,-.6,-.5,-.7,-.9,-.7,-.8, -0.9, -0.8, -0.9, -1.1, -1.2, -1.3]).reshape((-1, 1))    
#im2 = 0.2 + np.array([-0.8, -1.2, -0.9, -0.9, -1.1, -0.8, -1.2,
#                 -1.1, -1.3, -1.2, -1.3, -1.1, -1.2, -1.3]).reshape((-1, 1))
    
#im2 = 0.2 + np.array([-0.8, -1.2, -0.9, -0.9, -1.1, -0.8, -1.2,
#                 -0.7,-.75,-.75,-.8,-.9, -0.9,-.8,-.85, -0.8,-.85,-.8, -0.9, -1.1, -1.2, -1.3]).reshape((-1, 1))
    



   
     
noisy_sta_ind = np.arange(nsta)[obs_extra_sig[:,0]>0]

noisy_off_grid_inds = []
noisy_on_grid_inds = []
for i in range(nsta):
    if i in noisy_sta_ind:
        if i in off_grid_sta_ind:
            noisy_off_grid_inds.append(i)
        else:
            noisy_on_grid_inds.append(i)
noisy_off_grid_inds = np.array(noisy_off_grid_inds)
noisy_on_grid_inds = np.array(noisy_on_grid_inds)            
        



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


#%% Current Model
zeta = im2 - mu_im2
# Build the initial covariance matrix of the residuals and its inverse
sigma22 = phi_im2 * phi_im2.T
dist22 = np.abs(im2_dist.T - im2_dist)
# Weight the data with the extra sigma
J = np.ones_like(zeta)
J = phi_im2 / np.sqrt(phi_im2**2 + obs_extra_sig**2)
print(J)
corr_adj22 = J*J.T
np.fill_diagonal(corr_adj22, 1.0)

cov_im2im2_og1 = corr_adj22*lb_spatial_correlation(dist22,meas1_22,meas2_22)*sigma22
# Now add the extra uncertainty
#np.fill_diagonal(cov_im2im2_og1, cov_im2im2_og1.diagonal() + obs_extra_sig[:,0]**2)
cov_im2im2_inv = np.linalg.inv(cov_im2im2_og1)


# The event term and its variance (for the observation points)
cov_im2im2_dl = np.delete(np.delete(cov_im2im2_inv,elems,axis=0),elems,axis=1)
im2_dist_dl = np.delete(im2_dist,elems,axis=0)
z_12_dix = lb_spatial_correlation(np.zeros_like(im2_dist_dl),0.01*np.ones_like(im2_dist_dl),np.delete(im2_measure,elems,axis=0)) 
J_dl = J[dix].reshape([-1,1]) 
zeta_dl = zeta_dl       
bias_num = (z_12_dix*J_dl).T.dot(cov_im2im2_dl.dot(zeta_dl * (z_12_dix*J_dl)))
bias_denom = (z_12_dix*J_dl).T.dot(cov_im2im2_dl.dot(z_12_dix*J_dl))

sig_delta_B_im2_og = np.sqrt(1.0 / ((1.0 / tau_im2**2) + bias_denom))
delta_B_im2_og = sig_delta_B_im2_og**2 * bias_num



# The event term, its variance and the total within-event standard
# deviation (for the predictions)
sig_delta_B_im1_og = np.sqrt(1.0 / ((1.0 / tau_im1**2) + bias_denom))
delta_B_im1_og = sig_delta_B_im1_og**2 * bias_num

# Weight the data with the extra sigma
J = np.ones_like(zeta)
J = np.sqrt(phi_im2**2 +sig_delta_B_im2_og**2) / np.sqrt(phi_im2**2 +sig_delta_B_im2_og**2 + obs_extra_sig**2)

corr_adj22 = J*J.T
np.fill_diagonal(corr_adj22, 1.0)

# normalized residuals
debiased_zeta = J*(zeta - delta_B_im2_og)

# Distance matrices
dist11 = np.abs(im1_dist.T - im1_dist)
dist12 = np.abs(im1_dist - im2_dist.T)

             


# Covariance for im1 x im1 and im1 x im2
cov_im1im1_og = lb_spatial_correlation(dist11,0.01*np.ones_like(dist11),0.01*np.ones_like(dist11)) * ((phi_im1 * phi_im1.T) + (sig_delta_B_im1_og * sig_delta_B_im1_og.T))
cov_im1im2_og = (np.ones_like(phi_im1)*J.T)*(lb_spatial_correlation(dist12,0.01*np.ones_like(dist12),im2_measure.T*np.ones_like(im1_dist))*((phi_im1 * phi_im2.T) + (sig_delta_B_im1_og * sig_delta_B_im2_og.T)))
#cov_im1im1_phi = ga_spatial_correlation(dist11) * (phi_im1 * phi_im1.T)
#cov_im1im2_phi = Z_12*(np.ones_like(phi_im1)*J.T)*ga_spatial_correlation(dist12) * (phi_im1 * phi_im2.T)
#cov_im1im2_tau = Z_12*(np.ones_like(phi_im1)*J.T)*ga_spatial_correlation(dist12) * (sig_delta_B_im1_og  * sig_delta_B_im2_og .T)

# Build the updated (by bias) covariance matrix of the residuals and its inverse
cov_im2im2_og =  corr_adj22*lb_spatial_correlation(dist22,meas1_22,meas2_22) *((phi_im2 * phi_im2.T) + (sig_delta_B_im2_og * sig_delta_B_im2_og.T))
#cov_im2im2_phi = Z_22*ga_spatial_correlation(dist22) * (phi_im2 * phi_im2.T)
# Now add the extra uncertainty
#np.fill_diagonal(cov_im2im2_og, phi_im2**2 + sig_delta_B_im2_og**2 + obs_extra_sig**2)
#np.fill_diagonal(cov_im2im2_phi, phi_im2**2 + obs_extra_sig**2)
# Invert
cov_im2im2_inv = np.linalg.inv(cov_im2im2_og)
#cov_im2im2_phi_inv = np.linalg.inv(cov_im2im2_phi)

# Regression coefficient matrix
rcmatrix_og = cov_im1im2_og.dot(cov_im2im2_inv)
#rcmatrix_phi = cov_im1im2_phi.dot(cov_im2im2_phi_inv)

# Conditional normalized mean
mu_im1_im2_og = mu_im1 + rcmatrix_og.dot(debiased_zeta) + delta_B_im1_og

# Conditional normalized covariance
cov_im1im1_im2 = cov_im1im1_og - rcmatrix_og.dot(cov_im1im2_og.T)
#cov_im1im1_im2_phi = cov_im1im1_phi - rcmatrix_phi.dot(cov_im1im2_phi.T)

# Conditional standard deviation (sqrt of the diagonal of the
# conditional covariance matrix)
sig_im1im1_im2_diag_og = np.sqrt(np.clip(np.diag(cov_im1im1_im2), 0, np.inf))
#sig_im1im1_im2_phi_diag_og = np.sqrt(np.clip(np.diag(cov_im1im1_im2_phi), 0, np.inf))

#%%
fig = plt.figure(figsize=(10, 10))
txt = 'Figure 1. A synthetic, 1-dimensional example comparing the conditioned ground motion estimates of the proposed model and the W2018 model. (a) Unconditioned (𝝁_(𝒀_𝟏 )) and conditioned mean estimate (𝝁_(𝒀_𝟏 |𝒚_𝟐 )) of ground motion as a function of position (in log PGA in units of log g) for the W2018 and proposed updating methods. Also shown are the observed data and the biased mean estimated 𝝁_(𝒀_𝟏 )+𝝁_(𝑩_𝟏 |𝒚_𝟐 )for both updating methods. (b) Uncertainty in ground motion estimates in terms of the log standard deviation of PGA (in log g). Specifically, this subplot shows the unconditioned within and between-event standard deviations 𝝓_𝟏 and 𝝉_𝟏, respectively, as well as the standard deviation of the conditioned between-event uncertainty 𝝈_(𝑩_𝟏 |𝒚_𝟐 ) for the proposed and W2018 methods. (c) The updated total standard deviation of ground motions (in units of log g) 𝝈_(𝒀_𝟏 |𝒚_𝟐 ), as well as the standard deviation of the sum of the unconditioned within-event variance 𝝓_𝟏^𝟐 and the conditioned between-event variance 𝝈_(𝑩_𝟏 |𝒚_𝟐)^𝟐 for both the W2018 and proposed updating models.'
xlabel = r'\begin{center}X-axis\\*\textit{\small{' + txt + r'}}\end{center}'
plt.tight_layout(pad=3)
ax1 = fig.add_subplot(211)
ax1.plot(im1_dist, mu_im1,':k' ,label=r'GMPE $\mu_{Y_1}$', linewidth=2.1)

ax1.plot(im1_dist, mu_im1+delta_B_im1_og,'g--',alpha=1,linewidth=2.5, label=r'$\mu_{Y_1} +\mu_{B_1|y_2}$ (W2018)', zorder = 20)

ax1.plot(im1_dist, mu_im1+delta_B_im1,'brown',alpha=1,linewidth=2.5,label=r'$\mu_{Y_1} +\mu_{B_1|y_2}$ (Proposed)')
ax1.plot(im1_dist, mu_im1_im2_og,'--r' ,alpha=.7,label = r'$\mu_{Y_1|y_2}$ (W2018)',linewidth=1.6)
ax1.plot(im1_dist, mu_im1_im2, '-b',alpha=.4,label=r'$\mu_{Y_1|y_2}$ (Proposed)',linewidth=1.6)


#ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk',markersize=8, label=r'Off Grid $y_2$')
if noisy_on_grid_inds != []:
    ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$')
    ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$')
else:
    ax1.plot(im2_dist[on_grid_sta_ind], im2[on_grid_sta_ind], '*k', markersize = 10,label=r'On Grid $y_2$')

if noisy_off_grid_inds != []:
    ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red')


#ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
ax1.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.01, 1))    
ax1.set_ylabel('Ground Motions (log PGA)', fontsize=16)
ax1.tick_params(labelsize=16)
#ax1.set_ylim((-1.42, -0.76))
# ax1.set_ylim((-2.32, -0.65))



# Plot the standard deviation

# ax2 = fig.add_subplot(312)
# #ax2.plot(im1_dist, np.sqrt(phi_im1**2 + tau_im1**2), label='Prior total')
# ax2.plot(im1_dist, phi_im1,'m-.',alpha=1,label = r'Prior $\phi_1$',linewidth=2.1)

# ax2.plot(im1_dist, tau_im1,'k:', label=r'Prior $\tau_1$',linewidth=2.1)

# #ax2.plot(im1_dist, sig_im1im1_im2_diag_within,'m--',alpha=.8, label=r'Within $\sigma_{W_1|y_2}$',linewidth=2.5)
# #ax2.plot(im1_dist, sig_im1im1_im2_diag_between,'g-.',alpha=.6 ,label=r'Between $\sigma_{B_1|y_2}$',linewidth=2.5)
# ax2.plot(im1_dist, sig_delta_B_im1_og,'--g' ,alpha=1, label=r'$\sigma_{B_1|y_2}}$ (W2018)',linewidth=2.1 )
# ax2.plot(im1_dist, sig_delta_B_im1,'brown',alpha=1, label=r'$\sigma_{B_1|y_2}}$ (Proposed)',linewidth=2.1 )
# #ax2.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 2.1)
# #ax2.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (Proposed)',linewidth = 2.1)
# #ax2.plot(im1_dist, sig_im1im1_im2_diag_og,'--r' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=1.6)
# #ax2.plot(im1_dist, sig_im1im1_im2_diag,'-b',alpha=.4, label=r'$\sigma_{Y_1|y_2}$ (Proposed)',linewidth=1.6)

# #ax2.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=8)
# #if noisy_on_grid_inds != []:
# #    ax2.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy Data $\sigma_{y_2}$')
# #    ax2.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
# #else:
# #    ax2.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
# #
# #if noisy_off_grid_inds != []:
# #    ax2.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)
# #ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
# ax2.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.01, 1))#ax2.set_xlabel('Location (km)', fontsize=16)
# ax2.set_ylabel('Standard deviation (log PGA)', fontsize=16)
# ax2.tick_params(labelsize=16)
# ax2.set_ylim((-0.01, 0.75))

# Plot the standard deviation

ax3 = fig.add_subplot(212)
#ax2.plot(im1_dist, np.sqrt(phi_im1**2 + tau_im1**2), label='Prior total')
#ax3.plot(im1_dist, tau_im1,'k:', label=r'$\tau_1$',linewidth=2.1)

#ax2.plot(im1_dist, sig_im1im1_im2_diag_within,'m--',alpha=.8, label=r'Within $\sigma_{W_1|y_2}$',linewidth=2.5)
#ax2.plot(im1_dist, sig_im1im1_im2_diag_between,'g-.',alpha=.6 ,label=r'Between $\sigma_{B_1|y_2}$',linewidth=2.5)
#ax2.plot(im1_dist, sig_delta_B_im1, label='Conditional tau' )
ax3.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 2.5,zorder=20)
ax3.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (Proposed)',linewidth = 2.5)
ax3.plot(im1_dist, sig_im1im1_im2_diag_og,'--r' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=1.6)
ax3.plot(im1_dist, sig_im1im1_im2_diag,'-b',alpha=.4, label=r'$\sigma_{Y_1|y_2}$ (Proposed)',linewidth=1.6)

#x2.plot(im1_dist, phi_im1,label = 'Prior phi')
ax3.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=8)
if noisy_on_grid_inds != []:
    ax3.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy Data $\sigma_{y_2}$')
    ax3.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
else:
    ax3.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')

if noisy_off_grid_inds != []:
    ax3.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)
#ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
ax3.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.01, 1))
ax3.set_xlabel('Location (km)', fontsize=16)
ax3.set_ylabel('Standard deviation (log PGA)', fontsize=16)
ax3.set_ylim((-0.01, 0.75))
ax3.tick_params(labelsize=16)


fig.savefig('/Users/dengler/Documents/uncertainty paper figures/Figure1.jpg', format='jpg', dpi=1200, bbox_inches="tight")


fig = plt.figure(figsize=(10,6))
plt.plot([im1_dist[0],im1_dist[-1]],[0,0],'k')
plt.plot(im1_dist,delta_B_im1_og-delta_B_im1,linewidth=3)
plt.plot(im1_dist,mu_im1_im2_og-mu_im1_im2,linewidth=3)

fig = plt.figure(figsize=(10,6))

plt.plot(im1_dist, sig_im1im1_im2_diag_within,'r--',alpha=.7, label=r'$\sigma_{W_1|w_2}$',linewidth=1.6)
plt.plot(im1_dist, sig_im1im1_im2_diag_between,'b-.',alpha=.4 ,label=r'$\sqrt{Diag(cc^T)}σ_{Η|y_2}$',linewidth=2.1)
plt.plot(im1_dist, sig_im1im1_im2_diag,'-k',alpha=.4, label=r'$\sigma_{Y_1|y_2}$',linewidth=1.6)


#ax2.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 1.6)
#ax2.plot(im1_dist, sig_im1im1_im2_diag_og,'--r' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=2.8)

plt.plot(im1_dist, phi_im1,'g--',label = r'$\phi_1$',linewidth=2.1)
plt.plot(im1_dist, sig_delta_B_im1,'m',linestyle=':', label=r'$\sigma_{B_1|y_2}$',linewidth=2.1 )
plt.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$',linewidth = 1.5)

#ax2.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=8)
#if noisy_on_grid_inds != []:
#    ax2.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy Data $\sigma_{y_2}$')
#    ax2.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
#else:
#    ax2.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
#
#if noisy_off_grid_inds != []:
#    ax2.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)
#ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
plt.legend(loc='lower right',fontsize=10)
plt.xlabel('Location (km)', fontsize=16)
plt.ylabel('Standard deviation (log PGA)', fontsize=16)
plt.tick_params(labelsize=16)



#fig = plt.figure(figsize=(8, 5))
##ax1 = fig.add_subplot(211)
#plt.plot(im1_dist, mu_im1,':k' ,label=r'GMPE $\mu_{Y_1}$', linewidth=2.8)
#
##plt.plot(im1_dist, mu_im1+delta_B_im1_og,'g--',alpha=1,linewidth=1.6, label=r'$\mu_{Y_1} +\mu_{B_1|y_2}$ (W2018)')
#
#plt.plot(im1_dist, mu_im1+delta_B_im1,'brown',alpha=.7,linewidth=2.8,label=r'$\mu_{Y_1} +\mu_{B_1|y_2}$ ')
##plt.plot(im1_dist, mu_im1_im2_og,'--r' ,alpha=.7,label = r'$\mu_{Y_1|y_2}$ (W2018)',linewidth=2.8)
#plt.plot(im1_dist, mu_im1_im2, '-k',alpha=.4,label=r'$\mu_{Y_1|y_2}$',linewidth=2.8)
#
#
##ax1.plot(im1_dist, biased_mu, label='Actual biased mu')
#plt.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk',markersize=8, label=r'Off Grid $y_2$')
#if noisy_on_grid_inds != []:
#    plt.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$')
#    plt.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$')
#else:
#    plt.plot(im2_dist[on_grid_sta_ind], im2[on_grid_sta_ind], '*k', markersize = 10,label=r'On Grid $y_2$')
#
#if noisy_off_grid_inds != []:
#    plt.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red')
#
#
##ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
#plt.legend(loc='upper right',fontsize=14,bbox_to_anchor=(1.305, 1))    
#plt.ylabel('Ground Motions (log PGA)', fontsize=16)
#plt.tick_params(labelsize=16)
#plt.ylim((-1.42, -0.6))
#
##fig = plt.figure(figsize=(10,6))
##plt = fig.add_subplot(212)
#plt.plot(im1_dist,sig_im1im1_im2_diag/sig_im1im1_im2_diag_og)
#plt.plot(im1_dist,np.ones(len(im1_dist)),'k')
#
#fig = plt.figure(figsize=(8,5))
#plt.plot(im1_dist, phi_im1,'r--',alpha=.7,label = r'$\phi_1$',linewidth=2.8)
#
#
#plt.plot(im1_dist, sig_im1im1_im2_diag_within,'r-',alpha=.7, label=r'$\sigma_{W_1|w_2}$',linewidth=2.8)
#
#plt.plot(im1_dist, tau_im1,'b--',alpha=.7, label=r'$\tau_1$',linewidth=2.8)
#plt.plot(im1_dist, sig_delta_B_im1,'b:',linestyle=':',alpha=.7, label=r'$\sigma_{B_1|y_2}$',linewidth=2.8 )
#
#plt.plot(im1_dist, sig_im1im1_im2_diag_between,'b-',alpha=.7 ,label=r'$\sqrt{Diag(cc^T)}σ_{Η|y_2}$',linewidth=2.8)
##plt.plot(im1_dist, sig_im1im1_im2_diag,'-k',alpha=.4, label=r'$\sigma_{Y_1|y_2}$',linewidth=1.6)
#
#
##ax2.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 1.6)
##ax2.plot(im1_dist, sig_im1im1_im2_diag_og,'--r' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=2.8)
#
##plt.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$',linewidth = 1.5)
#
#plt.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=8)
#if noisy_on_grid_inds != []:
#    plt.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy Data $\sigma_{y_2}$')
#    plt.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
#else:
#    plt.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
#
#if noisy_off_grid_inds != []:
#    plt.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)
##ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
#plt.legend(loc='upper right',fontsize=14,bbox_to_anchor=(1.385, 1))
#plt.xlabel('Location (km)', fontsize=16)
#plt.ylabel('Standard deviation (log PGA)', fontsize=16)
#plt.tick_params(labelsize=16)

#fig=plt.figure(figsize=(10,6))
##ax1 = fig.add_subplot(211)
##ax1.plot(im1_dist,np.abs(mu_im1_im2)-np.abs(mu_im1)-(np.abs(mu_im1_im2_og)-np.abs(mu_im1)))
#plt.plot(im1_dist,(delta_B_im1)-((delta_B_im1_og)),label=r'(Bias (Prop))-(Bias (W2018))')
#plt.plot(im1_dist,mu_im1_im2-mu_im1_im2_og,label=r'($\mu_{Y_1|y_2}$ (Prop))-($\mu_{Y_1|y_2}$ (W2018))')
#plt.plot(im1_dist,np.zeros(len(im1_dist)),'k')
#plt.xlabel('Distance from Rupture (km)',fontsize=16)
#plt.ylabel('Proposed mean - W2018 mean',fontsize=16)
#plt.legend()
#
#sig_im1im1_im2_diag[sig_im1im1_im2_diag<1E-7]=.1
#sig_im1im1_im2_diag_og[sig_im1im1_im2_diag_og<1E-7]=.1
#
#fig=plt.figure(figsize=(10,6))



################
fig = plt.figure(figsize=(8.5, 8))
#ax1 = fig.add_subplot(211)
plt.plot(im1_dist, mu_im1,'k' ,label=r'GMM $\mu_{Y_1}$', linewidth=2.8)

plt.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk',markersize=11, label=r'Off Grid $y_2$')
if noisy_on_grid_inds != []:
    plt.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=13,label=r'Noisy $y_2$')
    plt.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 16,label=r'On Grid $y_2$')
else:
    plt.plot(im2_dist[on_grid_sta_ind], im2[on_grid_sta_ind], '*k', markersize = 18,label=r'On Grid $y_2$')

if noisy_off_grid_inds != []:
    plt.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red')

plt.plot(im1_dist, mu_im1+delta_B_im1_og,'g--',alpha=1,linewidth=2.8, label='Debiased Mean\n (W2018)')

plt.plot(im1_dist, mu_im1+delta_B_im1,'brown',linestyle='--',alpha=.7,linewidth=2.8,label='Debiased Mean')
plt.plot(im1_dist, mu_im1_im2_og,'g' ,alpha=1,label = r'$\mu_{Y_1|y_2}$ (W2018)',linewidth=2.8)
plt.plot(im1_dist, mu_im1_im2, 'brown',alpha=.7,label=r'$\mu_{Y_1|y_2}$',linewidth=2.8)


#ax1.plot(im1_dist, biased_mu, label='Actual biased mu')


#ax1.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
plt.legend(loc='upper left',fontsize=16,bbox_to_anchor=(1, 1),frameon=False)   
plt.xlabel('Location (km)',fontsize = 22) 
plt.ylabel('Ground Motions (log PGA)', fontsize=22)
plt.tick_params(labelsize=20)
plt.xlim([0,260])
plt.ylim((-1.42, -0.75))

#fig = plt.figure(figsize=(10,6))
#plt = fig.add_subplot(212)
#plt.plot(im1_dist,sig_im1im1_im2_diag/sig_im1im1_im2_diag_og)
#plt.plot(im1_dist,np.ones(len(im1_dist)),'k')



fig = plt.figure(figsize=(8.5,8))
#plt.plot(im1_dist, np.sqrt(tau_im1**2 + phi_im1**2),'k',label='Total',linewidth=2.8)
plt.plot(im1_dist, phi_im1,'r--',alpha=.8,label = r'$\phi_1$',linewidth=2.8)


#plt.plot(im1_dist, sig_im1im1_im2_diag_within,'r-',alpha=.7, label=r'$\sigma_{W_1|w_2}$',linewidth=2.8)

plt.plot(im1_dist, tau_im1,'b--',alpha=.8, label=r'$\tau_1$',linewidth=2.8)


plt.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=11)
if noisy_on_grid_inds != []:
    plt.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=13,label=r'Noisy Data $\sigma_{y_2}$')
    plt.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 16, label=r'On Grid $\sigma_{y_2}$')
else:
    plt.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 18, label=r'On Grid $\sigma_{y_2}$')

if noisy_off_grid_inds != []:
    plt.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)

plt.plot(im1_dist, sig_delta_B_im1_og,'g--',alpha=1, label=r'$\sigma_{B_1|y_2}$ (W2018)',linewidth=2.8 )
plt.plot(im1_dist, sig_delta_B_im1,'brown',linestyle='--',alpha=.7, label=r'$\sigma_{B_1|y_2}$',linewidth=2.8 )

#plt.plot(im1_dist, sig_im1im1_im2_diag_between,'b-',alpha=.7 ,label=r'$\sqrt{Diag(cc^T)}σ_{Η|y_2}$',linewidth=2.8)
#plt.plot(im1_dist, sig_im1im1_im2_diag,'-k',alpha=.4, label=r'$\sigma_{Y_1|y_2}$',linewidth=1.6)


#ax2.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 1.6)
plt.plot(im1_dist, sig_im1im1_im2_diag_og,'g' ,alpha=1, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=2.8)
plt.plot(im1_dist, sig_im1im1_im2_diag,'brown' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$',linewidth=2.8)

#plt.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$',linewidth = 1.5)

#ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
plt.legend(loc='upper left',fontsize=16,bbox_to_anchor=(1., 1),frameon=False)
plt.xlabel('Location (km)', fontsize=22)
plt.ylabel('Standard deviation (log PGA)', fontsize=22,labelpad =20)
plt.xlim([0,260])
plt.ylim([-.03,.75])
plt.tick_params(labelsize=20)


fig = plt.figure(figsize=(12,8))
#plt.plot(im1_dist, np.sqrt(tau_im1**2 + phi_im1**2),'k',label='Total',linewidth=2.8)
plt.plot(im1_dist, phi_im1,'r--',alpha=.8,label = r'$\phi_1$',linewidth=2.8)




plt.plot(im1_dist, tau_im1,'b--',alpha=.8, label=r'$\tau_1$',linewidth=2.8)


plt.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=11)
if noisy_on_grid_inds != []:
    plt.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=13,label=r'Noisy Data $\sigma_{y_2}$')
    plt.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 16, label=r'On Grid $\sigma_{y_2}$')
else:
    plt.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 18, label=r'On Grid $\sigma_{y_2}$')

if noisy_off_grid_inds != []:
    plt.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)

#plt.plot(im1_dist, sig_delta_B_im1_og,'g--',alpha=1, label=r'$\sigma_{B_1|y_2}$ (W2018)',linewidth=2.8 )
plt.plot(im1_dist, sig_delta_B_im1,'blue',linestyle=':',alpha=.7, label=r'$\sigma_{B_1|y_2}$',linewidth=2.8 )

plt.plot(im1_dist, sig_im1im1_im2_diag_within,'r-',alpha=.7, label=r'$\sigma_{W_1|w_2}$',linewidth=2.8)
plt.plot(im1_dist, sig_im1im1_im2_diag_between,'b-',alpha=.7 ,label=r'$\sqrt{Diag(cc^T)}σ_{Η|y_2}$',linewidth=2.8)
#plt.plot(im1_dist, sig_im1im1_im2_diag,'-k',alpha=.4, label=r'$\sigma_{Y_1|y_2}$',linewidth=1.6)


#ax2.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 1.6)
#plt.plot(im1_dist, sig_im1im1_im2_diag_og,'g' ,alpha=1, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=2.8)
#plt.plot(im1_dist, sig_im1im1_im2_diag,'brown' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=2.8)

#plt.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$',linewidth = 1.5)

#ax2.legend(numpoints=1, loc='lower left',bbox_to_anchor=(1.05, .5))
plt.legend(loc='upper left',fontsize=16,bbox_to_anchor=(1., 1),frameon=False)
plt.xlabel('Location (km)', fontsize=22)
plt.ylabel('Standard deviation (log PGA)', fontsize=22,labelpad =20)
plt.xlim([0,260])
plt.ylim([-.03,.75])
plt.tick_params(labelsize=20)