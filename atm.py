# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:38:00 2020

@author: Samuel
"""

import numpy as np
import pandas as pd


def getWind(h):
    
    wvec = np.array([0,0,0])
    
    return wvec

def getRho_from_table(atmdat, h):
    rho = np.interp(h, atmdat[0,:], atmdat[1,:])
    
    return rho

# def getRho_from_EarthGRAM(atmdat, h):
#     '''
#     gets density at given altitude from Earth-GRAM2016 formatted dataframe
#     '''
#     rho = np.interp(h, atmdat['Hgtkm'], atmdat['DensPert'])
    
#     return rho

def getMCAtmdat(filename, Nmc):
    '''
    loads atmosphere table for Nmc Monte Carlo trials. Written for a specially
    formatted output from Earth-GRAM2016.
    '''
    
    # Load one line at a time to find # of lines per MC run
    hprev = 1e8
    hcur = 1e7
    rMC = 0
    srows = []
    
    while hcur < hprev:
        atmRow = pd.read_csv(filename, delim_whitespace=True,
                              nrows=1, skiprows=srows)
        hprev = hcur
        hcur = atmRow['Hgtkm'][0]
        rMC += 1
        srows.append(rMC)
        
    # Now load requested number of profiles and return list of panda dfs
    Nloaded = 0
    datList = []
    while Nloaded < Nmc:
        atmdat = pd.read_csv(filename, delim_whitespace=True,
                             nrows = rMC-1,
                             skiprows=range(1, Nloaded*(rMC-1) + 1))
        # sort such that altitude is increasing
        atmdat.sort_values(by=['Hgtkm'], ascending=True, inplace=True)
        datList.append(atmdat)
        Nloaded += 1
        # print('loaded MC atm profile %d\n' %Nloaded)
    
    return datList

def getMCdens(filename, Nmc):
    '''
    loads atmosphere table for Nmc Monte Carlo trials. Written for a specially
    formatted output from Earth-GRAM2016.
    
    Similar to getMCAtmdat, but returns np array of just altitude, mean density,
    and each perturbed density.
    Ex:
        # ## get Nmc atmosphere profiles
        # Nmc = 1
        # i_trial = 0
        # densPert, densMean, h = getMCdens(filename, Nmc)
        # # at some point would be good to build this as a pandas df instead of np array
        # # rhoTable = np.array([h,densPert[:,i_trial]])
        # # load nominal density:
        # rhoTable = np.array([h,densPert[:,i_trial]])
        # params.atmdat = rhoTable
    '''
    
    # Load one line at a time to find # of lines per MC run
    hprev = 1e8
    hcur = 1e7
    rMC = 0
    srows = []
    
    while hcur < hprev:
        atmRow = pd.read_csv(filename, delim_whitespace=True,
                              nrows=1, skiprows=srows)
        hprev = hcur
        hcur = atmRow['Hgtkm'][0]
        rMC += 1
        srows.append(rMC)
        
    # Now load requested number of profiles and return numpy arrays
    Nloaded = 0
    densTot = np.empty([rMC-1, Nmc])
    densTot[:] = np.NaN
    
    while Nloaded < Nmc:
        atmdat = pd.read_csv(filename, delim_whitespace=True,
                             nrows = rMC-1,
                             skiprows=range(1, Nloaded*(rMC-1) + 1))
        # sort such that altitude is increasing
        atmdat.sort_values(by=['Hgtkm'], ascending=True, inplace=True)
        densTot[:,Nloaded] = atmdat['DensPert']
        Nloaded += 1
        print('loaded MC atm profile %d' %Nloaded)
        
    densMean = atmdat['DensMean']
    h = atmdat['Hgtkm']
    
    return densTot, densMean, h

def getMarsGRAMDensTable(filename, Nmc):
    '''
    loads atmosphere table for Nmc Monte Carlo trials. Written for a specially
    formatted output from Mars-GRAM2010.
    
    Similar to getMCAtmdat, but returns np array of just altitude, mean density,
    and each perturbed density.
    Ex:
        # ## get Nmc atmosphere profiles
        # Nmc = 1
        # i_trial = 0
        # densPert, densMean, h = getMCdens(filename, Nmc)
        # # at some point would be good to build this as a pandas df instead of np array
        # # rhoTable = np.array([h,densPert[:,i_trial]])
        # # load nominal density:
        # rhoTable = np.array([h,densPert[:,i_trial]])
        # params.atmdat = rhoTable
    '''
    
    # Load one line at a time to find # of lines per MC run
    hprev = -2
    hcur = -1
    rMC = 0
    srows = []
    
    while hcur > hprev:
        atmRow = pd.read_csv(filename, delim_whitespace=True,
                              nrows=1, skiprows=srows)
        hprev = hcur
        hcur = atmRow['HgtMOLA'][0]
        rMC += 1
        srows.append(rMC)
        
    # Now load requested number of profiles and return numpy arrays
    Nloaded = 0
    densTot = np.empty([rMC-1, Nmc])
    densTot[:] = np.NaN
    
    while Nloaded < Nmc:
        atmdat = pd.read_csv(filename, delim_whitespace=True,
                             nrows = rMC-1,
                             skiprows=range(1, Nloaded*(rMC-1) + 1))
        # sort such that altitude is increasing
        atmdat.sort_values(by=['HgtMOLA'], ascending=True, inplace=True)
        densTot[:,Nloaded] = atmdat['Denkgm3'] * atmdat['DensP']
        Nloaded += 1
        print('loaded MC atm profile %d' %Nloaded)
        
    densMean = atmdat['Denkgm3']
    h = atmdat['HgtMOLA']
    
    return densTot, densMean, h



