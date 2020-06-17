# -*- coding: utf-8 -*-
"""
Created on Sat May  2 11:53:51 2020

@author: Samuel
"""

import numpy as np
from scipy.integrate import solve_ivp

import ODEs

import sys


def simRun(params, tspan, events):
    # now do the setup and call dynamics
    y0 = np.block([params.x0, params.v0])
        
    sol = solve_ivp(lambda t, y: ODEs.dynamics(t,y,params), 
                    [tspan[0], tspan[-1]], y0.flatten(),
                    rtol=1e-9,atol=1e-9,
                    events=events) 
                    # normally set tols to 1e-12
    print(sol.message)
    
    if not sol.success:
        sys.exit('integration failed')
    
    if not sol.status:
        sys.exit('no termination event reached') 
    
    return sol