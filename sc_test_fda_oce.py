#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 18:17:56 2019

@author: roquet
"""
# Set the R_HOME environment variable to the R home directory
# To know the path of your R_HOME, execute in R the command RHOME
import os
os.environ['R_HOME'] = '//anaconda3/envs/r-conda/lib/R'

# Test if rpy2 is installed and print version
import rpy2
print("rpy2 version "+str(rpy2.__version__))

# change your working directory
path = "/Users/roquet/Documents/GitHub/fda_oce_python/"
os.chdir(path)

# load fda_oce_python functions
from fda_oce_python import *
import matplotlib.pyplot as plt

#%%
import scipy.io
data = scipy.io.loadmat('GLORYS_SO_2015-12.mat')
Pi=data['Pi']
Xi=data['Xi']

#%%
fdobj = bspl(Pi,Xi)
X_spli = data_from_fdobj(Pi,fdobj)

ii   = 599  #index of a profile
ndim = pca['ndim']

fig = plt.figure()
for kk in range(ndim): #Loop for each variable                        
    plt.subplot(1,ndim,kk+1)
    plt.plot(Xi[:,ii,kk],Pi,'.',X_spli[:,ii,kk],Pi)
    plt.ylim(plt.ylim()[::-1])
    plt.xlabel(pca['fdnames'][kk+2])
    if kk==0:
        plt.ylabel("Pressure")
    if kk==ndim-1:
        plt.legend(('raw','spline'))
    plt.suptitle("Example of a TS profile (#599)")

#%%
pca = fpca(fdobj)

#%%
pc = pc_from_fdobj(fdobj,pca)
pc_plot(pca, pc, sign=[-1,1])

#%%
eigenf_plot(pca, npc=1, sign=-1)


#%% Plots
fdobj_reco = fdobj_from_pc(pca,pc,2) # reconstruction with 2 modes
X_reco_1 = data_from_fdobj(Pi,fdobj_reco)

fdobj_reco = fdobj_from_pc(pca,pc,5) # reconstruction with 5 modes
X_reco_2 = data_from_fdobj(Pi,fdobj_reco)

ii   = 599  #index of a profile
ndim = pca['ndim']

fig = plt.figure()
for kk in range(ndim): #Loop for each variable                        
    plt.subplot(1,ndim,kk+1)
    plt.plot(Xi[:,ii,kk],Pi,'.',X_spli[:,ii,kk],Pi,X_reco_1[:,ii,kk],Pi,X_reco_2[:,ii,kk],Pi)
    plt.ylim(plt.ylim()[::-1])
    plt.xlabel(pca['fdnames'][kk+2])
    if kk==0:
        plt.ylabel("Pressure")
    if kk==ndim-1:
        plt.legend(('raw','spline','2 modes','5 modes'))
