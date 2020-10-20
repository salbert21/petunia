# -*- coding: utf-8 -*-
"""
UQ.py:
    contains functions related to KLE, PCE, and other UQ methods for petunia
    
Created on Tue Oct 20 14:25:03 2020

@author: Samuel Albert
"""


import numpy as np
from numpy import linalg as LA
from scipy.stats import norm


def getEproblem(filename, alpha, pfact):
    '''
    take in a filename and alpha, output the eigenvectors, eigenvalues, and d.
    INPUTS:
        filename - .npz file with density profiles
        alpha - accuracy parameter for truncating eigenvalues
        pfact - perturbation factor IN ADDITION to GRAM dispersions. 
            Set pfact = 1 to use GRAM data directly.
    OUTPUTS:
        evals - array of eigenvalues
        evecs - array of eigenvectors ordered same as evals array
        densSampMean - sample mean density profile, may not equal GRAM mean
        d - number of eigenvalues kept for density
        h - altitude profile array
    '''
    ## Load density data from the numpy binary file of EarthGRAM data
    densdata = np.load(filename)
    densTot = densdata['densTot']
    # densMean = densdata['arr_1'] # use sample mean instead of model mean!
    h = densdata['h']
    
    ## NOTE - we leave this out for now because we are going to assume density is
    #   Gaussian anyway, might as well calculate all evecs. In the future could
    #   turn this on to make data more truly Gaussian.
    # # TEMP - cuts off the density values near the surface ##
    # densTot = densTot[:200,:]
    # densMean = densMean[:200]
    # h = h[:200]
    # ##
    
    ## Get centered density about SAMPLE mean
    densSampMean = np.mean(densTot, axis=1)
    densCentered = densTot - densSampMean[:,None]
    
    densCentered = pfact * densCentered
    
    ## Get sample covariance matrix
    # this is equal to 1 / (N - 1) * densCentered @ densCentered.T (biased)
    covmat = np.cov(densCentered) # biased, by default, for sample estimate
    
    ## Eigenvalue decomposition of the sample covariance matrix, and sort desc.
    evals, evecs = LA.eig(covmat)
    idx = evals.argsort()[::-1] # descending order
    evals = evals[idx]
    evecs = evecs[:,idx]
    
    ## Find d, number of e-vals needed to get at least alpha accuracy
    rat = 0
    d = 0
    while rat < alpha:
        d += 1
        rat = np.sum(evals[0:d]) / np.sum(evals)
        
    return evals, evecs, densSampMean, d, h