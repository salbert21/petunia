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
    rho = np.interp(h, atmdat['alt'], atmdat['density'])
    
    return rho

def getRho_from_EarthGRAM(atmdat, h):
    '''
    gets density at given altitude from Earth-GRAM2016 formatted dataframe
    '''
    rho = np.interp(h, atmdat['Hgtkm'], atmdat['DensPert'])
    
    return rho

def getAtmdat(filename):
    atmdat = np.genfromtxt(filename, delimiter=',', names=True,
                           encoding='utf-8-sig')
    
    return atmdat

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
                             nrows = rMC-1, skiprows=range(1, Nloaded*rMC))
        # sort such that altitude is increasing
        atmdat.sort_values(by=['Hgtkm'], ascending=True, inplace=True)
        datList.append(atmdat)
        Nloaded += 1
    
    return datList


