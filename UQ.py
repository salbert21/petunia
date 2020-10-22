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
from scipy.special import factorial as fact


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



def getKLEdensfun(evals, evecs, mean, d, hvec):
    
    # # make sure h and mean are in ascending order for the interpolation
    # hvec = hvec[np.argsort(hvec)]
    # mean = mean[np.argsort(hvec)]
    
    # # first define function that evaluates mean at any h
    # def evalMean(h):
    #     return np.interp(h, hvec, mean)
        
    
    # generate d i.i.d. standard normal variables
    Ys = norm.rvs(size = d)
    rhovec = []
    
    # get density at each given altitude step
    for i in range(hvec.size):
        rhovec.append(mean[i] + np.sum(np.sqrt(evals[0:d]) * evecs[i,0:d] * Ys))
        
    rhovec = np.asarray(rhovec)
        
    # sort everything into ascending altitude order for interpolation
    rhovec = rhovec[np.argsort(hvec)]
    hvec = hvec[np.argsort(hvec)]
    
        
    # now define a function that simply interpolates rhovec for a given h
    def getKLErho(h):
        '''
        function for a single instance of the density KLE
        '''
        return np.interp(h, hvec, rhovec)

    return getKLErho, Ys


def nD_poly_array(d, p):
    '''
    Generates the multi-indices for all PC basis functions.
    Ported directly from MATLAB code written by Dr. A. Doostan.
    d: dimensions
    p: max total order
    '''
    
    MM = d-1
    n = 1
    
    t = np.empty(MM)
    t[:] = np.NaN
    
    PsiBasis = [[0] * d]
    
    for CurrentOrder in range(1,p+1):
        EndGenere = 0
        FirstThisOrder = 0
        
        while EndGenere == 0:
            n += 1
            # first list t for order CurrentOrder
            if FirstThisOrder == 0:
                for i in range(MM):
                    t[i] = i+1
                
                FirstThisOrder = 1
                
            else:
                # regular incrementation
                if t[MM-1] < (MM + CurrentOrder):
                    t[MM-1] += 1
                else:
                    j = MM
                    while t[j-1] == (j + CurrentOrder):
                        j -= 1
                    
                    t[j-1] += 1
                    for k in range(j+1, MM+1):
                        t[k-1] = t[j-1] + k - j
                        
            # direct translating t into PsiBasis(n)
            PsiBasis.append([]) # add blank new row to the list
            PsiBasis[n-1].append(t[0] - 1)
            # PsiBasis[(n-1),0] = t[0] - 1
            # print(n)
            for i in range (1,MM):
                PsiBasis[n-1].append(t[i] - t[i-1] - 1)
                # PsiBasis[n-1,i] = t[i] - t[i-1] - 1
            
            PsiBasis[n-1].append(d + CurrentOrder - t[MM-1] - 1)
            # PsiBasis[n-1,d-1] = d + CurrentOrder - t[MM-1] - 1
            
            # end of generation of order CurrentOrder
            if t[0] == (CurrentOrder + 1):
                EndGenere += 1
                
    PsiBasis = np.asarray(PsiBasis)
    
    PsiBasis = np.flip(PsiBasis, axis=1)
    
    return PsiBasis
            

def legendre_1d(p,x):
    '''
    evaluates 1D Legendre polynomial on [-1,1].
    Ported directly from MATLAB code written by Dr. A. Doostan.
    p is the max total order of the PC
    x is a d-dimensional array of eval points [-1,1]
    '''
    
    d = len(x)
    Leg = np.zeros([p+1, d])
    
    if p == 0:
        Leg[p,:] = np.ones(d)
    elif p == 1:
        Leg[p-1, :] = np.ones(d)
        Leg[p,:] = x
    else:
        Leg[0,:] = np.ones(d)
        Leg[1,:] = x
        for ord in range(2, p+1):
            Leg[ord,:] = (((2*ord)-1) * x * Leg[ord-1,:]\
                          - (ord-1) * Leg[ord-2,:])/ord
    
    # now normalize
    for i in range(p+1):
        Leg[i,:] = Leg[i,:] / np.sqrt(1 / (2*i + 1))
    
    return Leg

def hermite_1d(p,x):
    '''
    evaluates 1D Hermite polynomial on [-1,1].
    Ported directly from MATLAB code written by Dr. A. Doostan.
    p is the max total order of the PC
    x is a d-dimensional array of eval points [-1,1]
    '''
    d = len(x)
    Her = np.zeros([p+1, d])
    
    if p == 0:
        Her[p,:] = np.ones(d)
    elif p == 1:
        Her[p-1,:] = np.ones(d)
        Her[p,:] = x
    else:
        Her[0,:] = np.ones(d)
        Her[1,:] = x
        for ord in range(2,p+1):
            Her[ord,:] = x * Her[ord-1,:] - (ord-1) * Her[ord-2,:]
    
    # now normalize
    for i in range(p+1):
        Her[i,:] = Her[i,:] / np.sqrt(fact(i))
        
    return Her
    

def piset(xi, index_pc):
    '''
    evaluates a multi-D PCE basis at a given point (xi_1, ..., xi_d)
    Ported directly from MATLAB code written by Dr. A. Doostan.
    xi is the evaluation point
    index_pc is the multi-index for the PCE basis
    
    Note: right now this defaults to an all-Legendre basis. 
            can modify later to include Hermite, etc. for some/all d vars
    '''
    
    pc_xi = np.ones(index_pc.shape[0])
    
    p = np.sum(index_pc[index_pc.shape[0]-1,:])
    
    Legendre = legendre_1d(int(p), xi)
    # Hermite = hermite_1d(int(p), xi)
    
    for id in range(1,index_pc.shape[1]+1):
        nnz_index = np.where(index_pc[:,id-1] > 0)
        # nnz_index = nnz_index[0]
        if np.sum(nnz_index): # if at least one nonzero index
            index_pc_nnz = index_pc[nnz_index, id-1].flatten().astype(int)
            pc_xi[nnz_index] = pc_xi[nnz_index]\
                * Legendre[index_pc_nnz, id-1]
                # * Hermite[index_pc_nnz, id-1]
    
    return pc_xi