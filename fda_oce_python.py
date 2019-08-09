#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:22:16 2019

@author: roquet
"""

import numpy as np

#%% project profiles on a Bspline basis (use R-package fda)
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector,StrVector
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri
fda = importr('fda')
base = importr('base')


def create_basis(pmin,pmax,nbas=20):
    prange = FloatVector([pmin,pmax])
    breaks = FloatVector(prange[0]+(prange[1]-prange[0])*np.tan(np.linspace(0.,1.,nbas-2))/np.tan(1))
    breaks[0] = prange[0]
    breaks[-1] = prange[-1]
    basis  = fda.create_bspline_basis(rangeval = prange,nbasis = nbas,norder = 4,breaks = breaks)
    return basis


def bspl(Pi,Xi,nbas=20,fdn = ['Temperature','Salinity']):
    print("Converting the profiles into B-splines...")
    basis = create_basis(Pi[0],Pi[-1],nbas)
    with localconverter(numpy2ri.converter) as cv:
        Xi_R = cv.py2ro(Xi) # convert from *py*thon to *ro*bjects
    fdobj = fda.Data2fd(argvals = FloatVector(Pi),y = Xi_R,basisobj = basis,fdnames = StrVector(['Level','Station']+fdn))
    size = np.array(fdobj[0]).shape
    print("{0} B-splines computed for {1} variables.".format(size[1],size[2]))
    return fdobj


def get_coefs(fdobj):
    return np.array(fdobj[0])

def get_basis(fdobj):
    return fdobj[1]

def get_fdnames(fdobj):
    return list(fdobj[2])

def get_nbas(fdobj):
    return get_basis(fdobj)[3][0]

def get_rangeval(fdobj):
    return np.array(get_basis(fdobj)[2])

def compute_metric(fdobj):
    return fda.eval_penalty(get_basis(fdobj))


#%% compute PCA basis for a given set of profiles
def fpca(fdobj):
    
    coefs = get_coefs(fdobj)
    nbas = get_nbas(fdobj)
    nobs = coefs.shape[1]
    ndim = coefs.shape[2]
    metric = compute_metric(fdobj)
    
    # covariance matrix
    C = np.zeros((nobs,ndim*nbas))
    for kk in range(ndim):
        C[:,kk*nbas:(kk+1)*nbas] = np.squeeze(coefs[:,:,kk]).T

    Cm = np.mean(C,axis=0)
    Cc = C - Cm[np.newaxis,:]
    
    # inertia
    inertia = np.zeros(ndim)
    for kk in range(ndim):
        V = Cc[:,kk*nbas:(kk+1)*nbas].T @ Cc[:,kk*nbas:(kk+1)*nbas] @ metric / nobs
        inertia[kk] = np.trace( V )
    
    # metric M
    M = np.zeros((ndim*nbas,ndim*nbas))
    Mdeminv = M.copy()
    W = M.copy()
    aux = np.diag(np.ones(nbas))
    for kk in range(ndim):    
        M[kk*nbas:(kk+1)*nbas,kk*nbas:(kk+1)*nbas] = aux / inertia[kk]
        Mdeminv[kk*nbas:(kk+1)*nbas,kk*nbas:(kk+1)*nbas] = aux * np.sqrt(inertia[kk])
        W[kk*nbas:(kk+1)*nbas,kk*nbas:(kk+1)*nbas] = metric
    Mdem = np.sqrt(M);
    
    # metric W
    W = (W+W.T)/2.
    from numpy.linalg import cholesky, inv, eig
    Wdem = cholesky(W).T
    Wdeminv = inv(Wdem)
    
    # Cov matrix and eigen values/vectors
    V = Mdem @ Wdem @ Cc.T @ Cc @ Wdem.T @ Mdem /nobs
    pca_values, pca_vectors = eig(V)
    idx = pca_values.argsort()[::-1]
    pca_values = pca_values[idx]
    pca_vectors = pca_vectors[:,idx]
    
    # complete pca structure
    pca = {}
    
    pca['nbas'] = nbas
    pca['ndim'] = ndim
    
    pca['basis'] = get_basis(fdobj)
    pca['fdobj_ref'] = fdobj
    pca['fdnames'] = get_fdnames(fdobj)
    
    pca['Cm'] = Cm
    pca['inertia'] = inertia
    pca['W'] = W
    pca['M'] = M

    pca['values'] = pca_values
    pca['pval'] = 100 * pca_values / np.sum(pca_values)
    pca['vecnotWM'] = pca_vectors
    pca['vectors'] = Mdeminv @ Wdeminv @ pca_vectors
    pca['axes'] = pca['vectors'] * np.sqrt(pca['values'])[:,np.newaxis].T

    ## Comments
    Verif1 = pca['vectors'][:,0].T @ W @ M @ pca['vectors'][:,0] - 1.
    Verif2 = pca['vectors'][:,0].T @ W @ M @ pca['vectors'][:,1]
    if Verif1>1.e-10 or Verif2>1.e-10:
        print("Warning : The eigen values are not orthogonal.")
        print("Verif1 = "+str(Verif1))
        print("Verif2 = "+str(Verif2))
    
    return pca


#%% project B-spline representation of profiles on a pca basis
def pc_from_fdobj(fdobj,pca):
    
    coefs = get_coefs(fdobj)
    nbas  = pca['nbas']
    ndim  = pca['ndim']
    nobs  = coefs.shape[1]
    
    C = np.zeros((nobs,ndim*nbas))
    for kk in range(ndim):
        C[:,kk*nbas:(kk+1)*nbas] = np.squeeze(coefs[:,:,kk]).T
    Cc = C - pca['Cm'][np.newaxis,:]

    pc = Cc @ pca['W'].T @ pca['M'] @ pca['vectors']

    Verif3 = pca['values'][0] - pc[:,0].T @ pc[:,0] /nobs
    if Verif3>1.e-10:
        print("Warning : The eigen values are not equivalent to the PC norm.")
        print("Verif3 = "+str(pca['values'][0] - pc[:,0].T @ pc[:,0] /nobs))

    return pc


#%% compute B-spline coefficients from pc values (with trauncation)
def fdobj_from_pc(pca,pc,ntrunc=0):

    nbas  = pca['nbas']
    ndim  = pca['ndim']
    nobs  = pc.shape[0]

    if ntrunc==0:
        ntrunc = nbas*ndim

    coef = np.zeros((nbas,nobs,ndim))
    for kk in range(ndim):
        coef[:,:,kk] = pca['Cm'][kk*nbas:(kk+1)*nbas,np.newaxis] + \
            pca['vectors'][kk*nbas:(kk+1)*nbas,:ntrunc] @ pc[:,:ntrunc].T
        
    with localconverter(numpy2ri.converter) as cv:
        coef_R = cv.py2ro(coef) # convert from *py*thon to *ro*bjects
    fdobj = fda.fd(coef_R,pca['basis'],StrVector(pca['fdnames']))
    
    return fdobj


#%% compute data using B-spline coefficients
def data_from_fdobj(depth,fdobj):
    return np.array(fda.eval_fd(FloatVector(depth),fdobj))

    
#%% scatter plot of PC[n1]:PC[n2]
def pc_plot(pca, pc , npc=[1,2], sign=[1,1]):
    import matplotlib.pyplot as plt
    if len(npc)==2:
        fig, ax = plt.subplots()
        plt.scatter(sign[0]*pc[:,npc[0]-1],sign[1]*pc[:,npc[1]-1])
        plt.xlabel("PC{0} ({1:5.2f} %)".format(npc[0],pca['pval'][npc[0]]))
        plt.ylabel("PC{0} ({1:5.2f} %)".format(npc[1],pca['pval'][npc[1]]))
        plt.grid()


#%% Plot eigenfunctions 

def eigenf_plot(pca, npc=1, sign=1):
    import matplotlib.pyplot as plt
    nbas  = pca['nbas']
    ndim  = pca['ndim']
    basis = pca['basis']
    fdnames = StrVector(pca['fdnames'])
    Cm    = pca['Cm']
    prange = get_rangeval(pca['fdobj_ref'])
    depth  = np.arange(np.ceil(prange[0]),np.floor(prange[1]))

    fig = plt.figure()
    for kk in range(ndim):
        
        nn = npc - 1
        pb = np.round(100*np.sum(pca['vecnotWM'][kk*nbas:(kk+1)*nbas,nn]**2))    # Percentage of the bloc
        
        Cmean = Cm[kk*nbas:(kk+1)*nbas]
        Cplus = sign*pca['axes'][kk*nbas:(kk+1)*nbas,nn]
        
        fdobj_vp = fda.fd(FloatVector(Cmean + Cplus),basis,fdnames)
        fdobj_vm = fda.fd(FloatVector(Cmean - Cplus),basis,fdnames)
        fdobj_v0 = fda.fd(FloatVector(Cmean),basis,fdnames)
        
        v0 = data_from_fdobj(depth,fdobj_v0)
        vp = data_from_fdobj(depth,fdobj_vp)
        vm = data_from_fdobj(depth,fdobj_vm)
        
        plt.subplot(1,ndim,kk+1)
        plt.plot(v0, depth, vp, depth, vm, depth)
        plt.ylim(plt.ylim()[::-1])
        plt.xlabel("{0} ({1:4.1f} %)".format(pca['fdnames'][kk+2],pb))
        if kk==0:
            plt.ylabel("Pressure")

    plt.suptitle("PC{0} ({1:5.2f} %)".format(npc,pca['pval'][nn]))



