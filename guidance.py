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
from sim import Params
import planetaryconstants as constants

import ODEs
from conversions import getApses, getApsesSphPR


def getApoapsisError(ra, rp, params):
    '''
    For FNPAG Mode 1. Returns current error f(ts) = ra - raStar, km
    '''
    
    # if orbit is hyperbolic, make ra positive so error is well-behaved
    if ra < 0:
        print('ra < 0 in getApoapsisError function')
        return abs(ra) - params.raStar
    else:
        return ra - params.raStar
    
    
def getDVCost(ra, rp, params):
    '''
    For FNPAG Mode 2. Returns current total DV cost, km/s
    '''
    
    DV1 = np.sqrt(2 * params.p.mu)\
        * abs((np.sqrt(1/ra - 1/(ra + params.rpStar))\
               - np.sqrt(1/ra - 1 / (ra + rp))))
            
    DV2 = np.sqrt(2 * params.p.mu)\
        * abs((np.sqrt(1/params.rpStar - 1 / (params.raStar + params.rpStar))\
               - np.sqrt(1/params.rpStar - 1 / (ra + params.rpStar))))
        
    return DV1 + DV2


def getError(xx0vec, t, ts, sig0, sigd, phase, errFun, params):
    '''
    calls appropriate dynamics function and error function for FNPAG
    '''
    
    # make sure switching time has not already passed
    if phase == 1:
        if ts <= t:
            ts = t
            
    if phase == 1:
        xxvec = dynFNPAGPhase1(xx0vec, t, ts, sig0, sigd, params)
    elif phase == 2:
        xxvec = dynFNPAGPhase2(xx0vec, t, sigd, params)
    
    ra, rp = getApses(xxvec[:3,-1], xxvec[3:,-1], params)
    
    return errFun(ra, rp, params)

def getErrorSph(xx0vec, t, ts, sig0, sigd, phase, errFun, params):
    '''
    calls appropriate dynamics function and error function for FNPAG
    '''
    
    # make sure switching time has not already passed
    if phase == 1:
        if ts <= t:
            ts = t
            
    if phase == 1:
        xxvec = dynFNPAGPhase1Sph(xx0vec, t, ts, sig0, sigd, params)
    elif phase == 2:
        xxvec = dynFNPAGPhase2Sph(xx0vec, t, sigd, params)
    
    ra, rp = getApsesSphPR(xxvec[:,-1], params)
    
    return errFun(ra, rp, params)
        









def updateFNPAG(xxvec, t, ts, sig0, sigd, phase, mode, params,
                sphericalEOMs = True):
    '''
    Guidance update for FNPAG, used in predictor-corrector.
    Computes ts for phase 1, sigd for phase 2.
    Root-finds via Brent's Method for Mode 1, minimizes via Golden Section for
        Mode 2. Note that in SciPy the Brent method for minimize_scalar is an
        augmented version of Golden Section.
    INPUTS:
        xxvec: current cartesian inertial state, km, km/s
        t: current time, s
        phase: phase, 1 or 2
        mode: mode, 1 or 2 currently supported
        params: various constants, includes target state info
    OUTPUTS:
        MODE 1: tsi: updated switching time, s
        MODE 2: sigd: final bank angle, rad
    '''
    
    # update error function with current state and time values
    if mode == 1:
        errFun = getApoapsisError
    elif mode == 2:
        errFun = getDVCost
    else:
        sys.exit('Mode not recognized for FNPAG')
    
    if phase == 1:
        if sphericalEOMs:
                    getErr = lambda ts: getErrorSph(xxvec, t, ts, sig0, sigd,
                                     phase, errFun, params)
        else:     
            getErr = lambda ts: getError(xxvec, t, ts, sig0, sigd,
                                         phase, errFun, params)
        
        if params.errtol1 > 0:
            # if error already below tolerance, don't update control parameter
            if mode == 1 and abs(getErr(ts)) < params.errtol1:
                return ts
        
        
        # ts1 = ts - 5
        # ts2 = ts + 5
        
        # # make sure signs of f(a) and f(b) are different for Brent's Method
        # while mode == 1 and np.sign(getErr(ts1)) == np.sign(getErr(ts2)):
        #     ts1 -= 1
        #     ts2 += 1
        #     print(ts1)
        
        # print(getErr(ts1))
        # print(getErr(ts2))
        
        if mode == 1:
            # Brent's Method:
            tsi, res = brentq(getErr, params.ts1, params.ts2, full_output = True)
            converged = res.converged
            
            # # scipy automatic root-finding method:
            # sol = root_scalar(getErr, x0 = ts, x1 = ts+1)
            # tsi = sol.root
            # converged = sol.converged
            
        elif mode == 2:
            # Golden Section Method:
            res = minimize_scalar(getErr, bracket = (params.ts1, params.ts2),
                                  method = 'Brent')
            
            # # scipy minimize method:
            # res = minimize(getErr, x0 = ts)
            
            tsi = res.x
            converged = res.success
            
        else:
            sys.exit('Mode not recognized for FNPAG')
            
        if not converged:
            sys.exit('Failed to converge during guidance update'\
                     ' (phase {0:d}, mode {1:d})'.format(phase, mode))
        return tsi
    
    elif phase == 2:
        if sphericalEOMs:
                    getErr = lambda sigd: getErrorSph(xxvec, t, ts, sig0, sigd,
                                       phase, errFun, params)
        else:
            getErr = lambda sigd: getError(xxvec, t, ts, sig0, sigd,
                                           phase, errFun, params)
        
        if params.errtol2 > 0:
            # if error already below tolerance, don't update control parameter
            if mode == 1 and abs(getErr(sigd)) < params.errtol2:
                return sigd
        
        # sig1 = sigd - 1
        # sig2 = min(sigd + 1, 180)
        
        if mode == 1:
            # Brent's Method:
            sigdi, res = brentq(getErr, params.sig1, params.sig2, full_output = True)
            converged = res.converged
            
            # # scipy automatic root-finding method:
            # sol = root_scalar(getErr, x0 = sigd, x1 = sigd+0.02)
            # sigdi = sol.root
            # converged = sol.converged
            
        elif mode == 2:
            # Golden Section method:
            res = minimize_scalar(getErr, bracket = (params.sig1, params.sig2),
                                  method = 'Brent')
            
            # # scipy minimize method:
            # res = minimize(getErr, x0 = sigd)
            
            sigdi = res.x
            converged = res.success
            
        else:
            sys.exit('Mode not recognized for FNPAG')
            
        if not converged:
            sys.exit('Failed to converge during guidance update'\
                     ' (phase {0:d}, mode {1:d})'.format(phase, mode))
        return sigdi
    
    
    
    
    
    
    
    
    
    
def dynFNPAGPhase1Sph(xx0vec, t, ts, sig0, sigd, params, returnTime = False):
    '''
    Dynamics function for Phase 1 of FNPAG.
    Propagages from given state to atm exit. Switches bank angle from sig0 to
        sigd instantaneously at time ts.
    INPUTS:
        xx0vec: initial spherical planet-relative state, m, m/s, rad
        t: initial (current) time, s
        ts: switching time, s
        sig0: initial bank angle, deg
        sigd: final bank angle, deg
        params: various constants
    OUTPUTS:
        xxvec: cartesian inertial state vector over time, km, km/s
    '''
    
    event1 = lambda t, y: ODEs.above_max_alt_sph(t, y, params)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt_sph(t, y, params)
    event2.terminal = True
    
    # propagate lift-up until ts
    tspan = (t, ts)
    sol1 = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMs(t, y,
                                                          np.radians(sig0),
                                                          params),
                     tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                     events = (event1, event2))
    if sol1.status != 0:
        sys.exit('switching time never reached during phase 1 of FNPAG')
    
    # propagate lift-down until atm exit or surface impact
    tspan = (ts, params.tf)
    sol2 = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMs(t, y,
                                                          np.radians(sigd),
                                                          params),
                     tspan, sol1.y[:,-1], rtol = params.rtol, atol = params.atol,
                     events = (event1, event2))
    
    if sol2.status != 1:
        sys.exit('never reached a terminal condition during phase 1 of FNPAG')
    
    xxvec = np.hstack([sol1.y, sol2.y])
    tvec = np.hstack([sol1.t, sol2.t])
    
    if returnTime:
        return xxvec, tvec
    else:
        return xxvec
    
def dynFNPAGPhase2Sph(xx0vec, t, sigd, params, returnTime = False):
    '''
    Dynamics function for Phase 2 of FNPAG.
    Propagates from given state to atm exit. Assumes constant bank angle.
    INPUTS:
        xx0vec: initial spherical planet-relative state, km, km/s
        t: initial (current) time, s
        sigd: final bank angle, deg
    OUTPUTS:
        xxvec: cartesian inertial state vector over time, km, km/s
    '''
    
    event1 = lambda t, y: ODEs.above_max_alt_sph(t, y, params)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt_sph(t, y, params)
    event2.terminal = True
    
    # propagate until tf or terminal condition
    tspan = (t, params.tf)
    sol = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMs(t, y,
                                                         np.radians(sigd),
                                                         params),
                    tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                    events = (event1, event2))
    if sol.status != 1:
        sys.exit('never reached a terminal condition during phase 2 of FNPAG')
    
    if returnTime:
        return sol.y, sol.t
    else:
        return sol.y
    
    
    
    


def dynFNPAGPhase1(xx0vec, t, ts, sig0, sigd, params, returnTime = False):
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
    sol1 = solve_ivp(lambda t, y: ODEs.dynamics(t, y, sig0, params),
                     tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                     events = (event1, event2))
    if sol1.status != 0:
        sys.exit('switching time never reached during phase 1 of FNPAG')
    
    # propagate lift-down until atm exit or surface impact
    tspan = (ts, params.tf)
    sol2 = solve_ivp(lambda t, y: ODEs.dynamics(t, y, sigd, params),
                     tspan, sol1.y[:,-1], rtol = params.rtol, atol = params.atol,
                     events = (event1, event2))
    
    if sol2.status != 1:
        sys.exit('never reached a terminal condition during phase 1 of FNPAG')
    
    xxvec = np.hstack([sol1.y, sol2.y])
    tvec = np.hstack([sol1.t, sol2.t])
    
    if returnTime:
        return xxvec, tvec
    else:
        return xxvec

def dynFNPAGPhase2(xx0vec, t, sigd, params, returnTime = False):
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
    sol = solve_ivp(lambda t, y: ODEs.dynamics(t, y, sigd, params),
                    tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                    events = (event1, event2))
    if sol.status != 1:
        sys.exit('never reached a terminal condition during phase 2 of FNPAG')
    
    if returnTime:
        return sol.y, sol.t
    else:
        return sol.y
    
    