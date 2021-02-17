# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:24:45 2021

@author: Samuel Albert
"""

from sim import Params
import planetaryconstants as constants
from atm import getRho_from_table

import numpy as np
import scipy.interpolate as interp
import time


### CREATE params INPUT CLASS FOR NOMINAL VALUES
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
params.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
params.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])

scipyfun = interp.interp1d(params.atmdat[0,:], params.atmdat[1,:])

def myinterp(h, params):
    '''
    h is altitude in km. assumes atm data sorted in ascending order!
    also, will not work if h is outside of data range!
    '''
    
    # get index of first element in altitude array great than h
    ind = np.argmax(params.atmdat[0,:] > h)
    x0 = params.atmdat[0, ind-1]
    x1 = params.atmdat[0, ind]
    y0 = params.atmdat[1, ind-1]
    y1 = params.atmdat[1, ind]
    
    rho = (y0 * (x1 - h) + y1 * (h - x0)) / (x1 - x0)
    
    return rho

h = 96.15498715
tic1 = time.time()
rho1 = getRho_from_table(params.atmdat, h)
toc1 = time.time()

tic2 = time.time()
rho2 = myinterp(h, params)
toc2 = time.time()

tic3 = time.time()
rho3 = scipyfun(h)
toc3 = time.time()

print(rho1)
print(rho2)
print(rho3)

# NOTE: just do timeit _____ for each rho command, in the IPython prompts

# print('numpy took {:.4e} s'.format(toc1-tic1))
# print('custom took {:.4e} s'.format(toc2-tic2))
# print('scipy took ')


