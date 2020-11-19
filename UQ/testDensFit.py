# -*- coding: utf-8 -*-
"""
testDensFit.py:
    check Gaussianity of density from MarsGRAM
Created on Wed Nov 18 18:03:49 2020

@author: Samuel Albert
"""

import numpy as np
# from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
from scipy.stats import norm
# import scipy.optimize as op
import time
from functools import partial
from statsmodels.graphics.gofplots import qqplot

# Define narrower bandwidth for kernel density estimation using function
#   straight from scipy.stats reference
def my_kde_bandwidth(obj, fac=1./5):
    """We use Scott's Rule, multiplied by a constant factor."""
    return np.power(obj.n, -1./(obj.d+4)) * fac

tic = time.time()

plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 15})

# =============================================================================
# Load density data
# =============================================================================
filename = '../data/Mars_0.1_50000.npz'
densdata = np.load(filename)
densTot = densdata['densTot']
h = densdata['h']
densMean = densdata['densMean']
densCentered = densTot - densMean[:,None]


# =============================================================================
# Loop through altitudes of interest
# =============================================================================
# for i in range(10):
# for i in np.linspace(-1,-500,10).astype(int):
for i in range(0,1430,130):
    if i == 1300:
        i = 1299

    # =========================================================================
    # kde vs Gaussian vs histogram
    # =========================================================================
    kde = stats.gaussian_kde(densCentered[i,:],
                          bw_method=partial(my_kde_bandwidth, fac=0.5))
    
    x_eval = np.linspace(densCentered[i,:].min(), densCentered[i,:].max(), 1000)
    
    params = stats.describe(densCentered[i,:])
    dmean = params[2]
    dvar = params[3]
    dstd = dvar**(1/2)
    
    rv = norm(loc=dmean, scale=dstd)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_eval, kde(x_eval), label='kde pdf', linewidth=4)
    ax.plot(x_eval, rv.pdf(x_eval), label='Gaussian pdf', linewidth=4)
    ax.hist(densCentered[i,:], bins='auto', density=True, label='histogram')
    
    ax.legend()
    ax.set_xlabel('Centered density, kg/m^3')
    ax.set_ylabel('pdf value at this density value/bin')
    ax.grid()
    ax.set_title('{} km altitude'.format(h[i]))
    
    
    # =========================================================================
    # q-q plot
    # =========================================================================
    qqplot(densCentered[i,:], line='s')
    plt.grid()
    plt.title('{} km altitude'.format(h[i]))















