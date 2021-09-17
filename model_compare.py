#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:09:06 2020

@author: dengler, cworden
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, uniform
import random
from numpy.linalg import cholesky


# %%-------------------------------------------------------------------------
# Loth & Baker (2013) spatial correlation for PGA
# Input: a numpy array of distances in km
# Output: spatial correlation coefficients (from 0.0 to 1.0)

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

         
# %%-------------------------------------------------------------------------

# Notation:
# im1 = predictions
# im2 = observations
# phi = within-event standard deviation
# tau = between-event standard deviation
# Number of points in the X_1 (prediction) vector
ngrid = 261
sgrid = 200

# Flag to determine if the usual toy problem is used or simulate new data
rand_flag = 1

# Number of stations with data
if rand_flag == 0:
    nsta = 14 
else
    nsta = 22

# Determines the distance from rupture to ignore station contributions to bias, and flag deciding whether to do this or not
cutoff_dist = 100
cutoff_flag =0

# Locations (in km) of the points in X_1
im1_dist = np.linspace(0.0, 261.0, ngrid).reshape((-1, 1))

if rand_flag == 0:
    ## Original toy problem 
    ind_obs = np.array([10, 20, 30, 40, 44, 48, 52,
                 90, 98, 106, 114, 160, 168, 176], dtype=int)
else:    
    # Select random observation indices
    ind_obs = np.sort(np.array(random.sample(range(sgrid),nsta)))    

nsta = len(ind_obs)

# Locations of stations
im2_dist = im1_dist[ind_obs]


# Select off-grid points
if rand_flag == 0:
    off_grid_sta_ind = np.arange(nsta)[-3:]
    im2_dist[[off_grid_sta_ind]] += 0.5 
else: 
# (randomly selects a random number of stations and distributes
# them uniformly between the grid points)
    num_off_grid_sta = random.sample(range(1,int(nsta/2)+1),1)[0]
    off_grid_sta_ind = np.sort(np.array(random.sample(range(nsta),num_off_grid_sta)))
    im2_dist[off_grid_sta_ind] += uniform.rvs(size=(num_off_grid_sta,1))

#Indices of the stations which are on grid
on_grid_sta_ind = np.delete(np.arange(nsta),off_grid_sta_ind)

# "GMPE" estimates at the locations in the X_1 array; here we use
# a constant (or a linear ramp) for demonstration purposes
mu_im1 = np.full_like(im1_dist, -1.0)

# "GMPE" phi and tau at locations in X_1; (increasing with distance for
# demonstration purposes)
#phi_im1 = 0.4 + 0.001*im1_dist      #Linearly increasing phi
phi_im1 = np.full_like(im1_dist, 0.6)     # Constant phi
tau_im1 = 0.25 + 0.001*im1_dist      #Linearly increasing tau
#tau_im1 = np.full_like(im1_dist, 0.3)     # Constant tau

mu_im2 = mu_im1[ind_obs]
phi_im2 = phi_im1[ind_obs]
tau_im2 = tau_im1[ind_obs]

# Additional sigma for some of the observations
if rand_flag  == 0 :
    ###Toy problem extra noise
    obs_extra_sig = np.array([0.3, 0, 0, 0, 0, 0, 0,
                           0.5, 0.5, 0.5, 0.5, 0, 0, 0]).reshape(-1, 1)
else:
    obs_extra_sig = np.full_like(phi_im2,0.0) #No noise
    num_obs_extra_sig = random.sample(range(1,int(nsta/2)+1),1)[0]
    noisy_sta_ind = np.sort(np.array(random.sample(range(nsta),num_obs_extra_sig)))
    obs_extra_sig[noisy_sta_ind] = .5

# Intensity measure of data
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
between_res = norm.rvs(size=1)
biased_mu = mu_im1 + tau_im1*between_stand_norm_real
between_res = between_stand_norm_real*tau_im2
#np.random.seed(11)
within_res = cholesky(cov_dW2dW2).dot(norm.rvs(size=(nsta,1)))

# Observations
if rand_flag == 0:
    #Original toy problem data
    im2 =-.1+  np.array([-.8, -1.2, -0.9, -0.9, -1.1, -0.8, -1.2,
                 -.7, -0.9, -0.8, -0.9, -1.1, -1.2, -1.3]).reshape((-1, 1))
else:
    im2 = mu_im2 + between_res + within_res 
     
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

# Build the updated (by bias) covariance matrix of the residuals and its inverse
cov_im2im2_og =  corr_adj22*lb_spatial_correlation(dist22,meas1_22,meas2_22) *((phi_im2 * phi_im2.T) + (sig_delta_B_im2_og * sig_delta_B_im2_og.T))

# Invert
cov_im2im2_inv = np.linalg.inv(cov_im2im2_og)

# Regression coefficient matrix
rcmatrix_og = cov_im1im2_og.dot(cov_im2im2_inv)

# Conditional normalized mean
mu_im1_im2_og = mu_im1 + rcmatrix_og.dot(debiased_zeta) + delta_B_im1_og

# Conditional normalized covariance
cov_im1im1_im2 = cov_im1im1_og - rcmatrix_og.dot(cov_im1im2_og.T)

# Conditional standard deviation (sqrt of the diagonal of the
# conditional covariance matrix)
sig_im1im1_im2_diag_og = np.sqrt(np.clip(np.diag(cov_im1im1_im2), 0, np.inf))

#%%

#### Plot the three panel figure summarizing the differences between the proposed and current W2018 conditioning method 
# for the 1-dimensional example
fig = plt.figure(figsize=(10, 15))
plt.tight_layout(pad=3)

# Plot the mean ground motions
ax1 = fig.add_subplot(311)
ax1.plot(im1_dist, mu_im1,':k' ,label=r'GMPE $\mu_{Y_1}$', linewidth=2.1)
ax1.plot(im1_dist, mu_im1+delta_B_im1_og,'g--',alpha=1,linewidth=2.5, label=r'$\mu_{Y_1} +\mu_{B_1|y_2}$ (W2018)', zorder = 20)
ax1.plot(im1_dist, mu_im1+delta_B_im1,'brown',alpha=1,linewidth=2.1,label=r'$\mu_{Y_1} +\mu_{B_1|y_2}$ (Proposed)')
ax1.plot(im1_dist, mu_im1_im2_og,'--r' ,alpha=.7,label = r'$\mu_{Y_1|y_2}$ (W2018)',linewidth=1.6)
ax1.plot(im1_dist, mu_im1_im2, '-b',alpha=.4,label=r'$\mu_{Y_1|y_2}$ (Proposed)',linewidth=1.6)
ax1.plot(im2_dist[off_grid_sta_ind], im2[off_grid_sta_ind], 'sk',markersize=8, label=r'Off Grid $y_2$')
if noisy_on_grid_inds != []:
    ax1.plot(im2_dist[noisy_on_grid_inds],im2[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy $y_2$')
    ax1.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], im2[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10,label=r'On Grid $y_2$')
else:
    ax1.plot(im2_dist[on_grid_sta_ind], im2[on_grid_sta_ind], '*k', markersize = 10,label=r'On Grid $y_2$')

if noisy_off_grid_inds != []:
    ax1.plot(im2_dist[noisy_off_grid_inds],im2[noisy_off_grid_inds],'sk',color='red')
ax1.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.01, 1))    
ax1.set_ylabel('Ground Motions (log PGA)', fontsize=16)
ax1.tick_params(labelsize=16)
ax1.set_ylim((-1.42, -0.76))


# Plot the prior standard deviations and bias standard deviations
ax2 = fig.add_subplot(312)
ax2.plot(im1_dist, phi_im1,'m-.',alpha=1,label = r'Prior $\phi_1$',linewidth=2.1)
ax2.plot(im1_dist, tau_im1,'k:', label=r'Prior $\tau_1$',linewidth=2.1)
ax2.plot(im1_dist, sig_delta_B_im1_og,'--g' ,alpha=1, label=r'$\sigma_{B_1|y_2}}$ (W2018)',linewidth=2.1 )
ax2.plot(im1_dist, sig_delta_B_im1,'brown',alpha=1, label=r'$\sigma_{B_1|y_2}}$ (Proposed)',linewidth=2.1 )
ax2.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.01, 1))#ax2.set_xlabel('Location (km)', fontsize=16)
ax2.set_ylabel('Standard deviation (log PGA)', fontsize=16)
ax2.tick_params(labelsize=16)
ax2.set_ylim((-0.01, 0.75))

# Plot the updated standard deviation
ax3 = fig.add_subplot(313)
ax3.plot(im1_dist, np.sqrt(sig_delta_B_im1_og**2 +phi_im1**2),'g--',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (W2018)',linewidth = 2.5,zorder=20)
ax3.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$ (Proposed)',linewidth = 2.1)
ax3.plot(im1_dist, sig_im1im1_im2_diag_og,'--r' ,alpha=.7, label=r'$\sigma_{Y_1|y_2}$ (W2018)',linewidth=1.6)
ax3.plot(im1_dist, sig_im1im1_im2_diag,'-b',alpha=.4, label=r'$\sigma_{Y_1|y_2}$ (Proposed)',linewidth=1.6)
ax3.plot(im2_dist[off_grid_sta_ind], obs_extra_sig[off_grid_sta_ind], 'sk', label=r'Off Grid $\sigma_{y_2}$',markersize=8)
if noisy_on_grid_inds != []:
    ax3.plot(im2_dist[noisy_on_grid_inds],obs_extra_sig[noisy_on_grid_inds],'^k',markersize=10,label=r'Noisy Data $\sigma_{y_2}$')
    ax3.plot(im2_dist[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], obs_extra_sig[np.delete(on_grid_sta_ind,noisy_on_grid_inds)], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
else:
    ax3.plot(im2_dist[on_grid_sta_ind], obs_extra_sig[on_grid_sta_ind,], '*k', markersize = 10, label=r'On Grid $\sigma_{y_2}$')
if noisy_off_grid_inds != []:
    ax3.plot(im2_dist[noisy_off_grid_inds],obs_extra_sig[noisy_off_grid_inds],'sk',color='red',markersize=8)
ax3.legend(loc='upper left',fontsize=12, bbox_to_anchor=(1.01, 1))
ax3.set_xlabel('Location (km)', fontsize=16)
ax3.set_ylabel('Standard deviation (log PGA)', fontsize=16)
ax3.set_ylim((-0.01, 0.75))
ax3.tick_params(labelsize=16)


#fig.savefig('model_comp.jpg', format='jpg', dpi=1200, bbox_inches="tight")


## Make a figure showing the separated prior and updated residual standard deviations for the 1-dimensional case
fig = plt.figure(figsize=(10,6))
plt.plot(im1_dist, sig_im1im1_im2_diag_within,'r--',alpha=.7, label=r'$\sigma_{W_1|w_2}$',linewidth=1.6)
plt.plot(im1_dist, sig_im1im1_im2_diag_between,'b-.',alpha=.4 ,label=r'$\sqrt{Diag(cc^T)}σ_{Η|y_2}$',linewidth=2.1)
plt.plot(im1_dist, sig_im1im1_im2_diag,'-k',alpha=.4, label=r'$\sigma_{Y_1|y_2}$',linewidth=1.6)
plt.plot(im1_dist, phi_im1,'g--',label = r'$\phi_1$',linewidth=2.1)
plt.plot(im1_dist, sig_delta_B_im1,'m',linestyle=':', label=r'$\sigma_{B_1|y_2}$',linewidth=2.1 )
plt.plot(im1_dist, np.sqrt(sig_delta_B_im1**2 +phi_im1**2),'brown',alpha=1,label=r'$\sqrt{\phi_1^2 + \sigma^2_{B_1|y_2}}$',linewidth = 1.5)
plt.legend(loc='lower right',fontsize=10)
plt.xlabel('Location (km)', fontsize=16)
plt.ylabel('Standard deviation of log PGA', fontsize=16)
plt.tick_params(labelsize=16)

