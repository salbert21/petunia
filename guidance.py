# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:01:58 2021

guidance.py:
    guidance-related functions for Petunia, including FNPEG and FNPAG algs

@author: Samuel Albert
"""


import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar, root_scalar, minimize
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

import ODEs


def dynFNPAGPhase1(xx0vec, t, ts, sig0, sigd, params):
    '''
    Dynamics function for Phase 1 of FNPAG.
    Propagages from given state to atm exit. Switches bank angle from sig0 to
        sigd instantaneously at time ts.
    INPUTS:
        xx0vec: initial cartesian inertial state, km, km/s
        t: initial (current) time, s
        ts: switching time, s
        sig0: initial bank angle, deg
        sigd: final bank angle, deg
        params: various constants
    OUTPUTS:
        xxvec: cartesian inertial state vector over time, km, km/s
    '''
    
    event1 = lambda t, y: ODEs.above_max_alt(t, y, params)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt(t, y, params)
    event2.terminal = True
    
    # propagate lift-up until ts
    tspan = (t, ts)
    params.bank = sig0
    sol1 = solve_ivp(lambda t, y: ODEs.dynamics(t, y, params),
                     tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                     events = (event1, event2))
    if sol1.status != 0:
        sys.exit('switching time never reached during phase 1 of FNPAG')
    
    # propagate lift-down until atm exit or surface impact
    tspan = (ts, params.tf)
    params.bank = sigd
    sol2 = solve_ivp(lambda t, y: ODEs.dynamics(t, y, params),
                     tspan, sol1.y[:,-1], rtol = params.rtol, atol = params.atol,
                     events = (event1, event2))
    
    if sol2.status != 1:
        sys.exit('never reached a terminal condition during phase 1 of FNPAG')
    
    xxvec = np.hstack([sol1.y, sol2.y])
    
    return xxvec

def dynFNPAGPhase2(xx0vec, t, sigd, params):
    '''
    Dynamics function for Phase 2 of FNPAG.
    Propagates from given state to atm exit. Assumes constant bank angle.
    INPUTS:
        xx0vec: initial cartesian inertial state, km, km/s
        t: initial (current) time, s
        sigd: final bank angle, deg
    OUTPUTS:
        xxvec: cartesian inertial state vector over time, km, km/s
    '''
    
    event1 = lambda t, y: ODEs.above_max_alt(t, y, params)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt(t, y, params)
    event2.terminal = True
    
    # propagate until tf or terminal condition
    tspan = (t, params.tf)
    params.bank = sigd
    sol = solve_ivp(lambda t, y: ODEs.dynamics(t, y, params),
                    tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                    events = (event1, event2))
    if sol.status != 1:
        sys.exit('never reached a terminal condition during phase 2 of FNPAG')
    
    return sol.y
    
    