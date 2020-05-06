# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:52:56 2020

@author: Samuel
"""

import numpy as np
from frames import DCMs
from atm import getRho_from_table, getWind

def dynamics(t,yy,params):
    '''
    ODEs for full dynamics acting on the vehicle
    assume:
        - just one vehicle
        - yy is [position; velocity] in INERTIAL frame
    '''
       
    
    dydt = np.empty(6)
    dydt[:] = np.NaN
    
    # extract inertial state components
    x = yy[0]
    y = yy[1]
    z = yy[2]
    dx = yy[3]
    dy = yy[4]
    dz = yy[5]
    
    xvec_N = np.array([x, y, z])
    vvec_N = np.array([dx, dy, dz])
    
    r = np.linalg.norm(xvec_N)
    mu_r3 = params.mu / r**3
    
    ## Get gravitational force
    Fgvec_N = - mu_r3 * params.m * xvec_N
    
    ## Get aerodynamics forces
    SN = DCMs.getSN(t, params)
    NS = DCMs.getNS(t, params)
    
    OMvec = np.array([0, 0, -params.p.om])
    
    vvec_S = SN @ (vvec_N + np.cross(OMvec, xvec_N))
    
    h = r - params.p.rad
    
    if params.dMode == 'fun':
        rho = params.dFun(h)
    else:
        rho = getRho_from_table(params.atmdat, h)
    wvec_S = getWind(h)
    
    vInfvec_S = vvec_S + wvec_S
    vInf = np.linalg.norm(vInfvec_S)
    
    vInfvec_N = NS @ vInfvec_S
    
    Lmag = 1/2 * rho * (vInf*1e3)**2 * params.CL * params.A / 1e3
    Dmag = 1/2 * rho * (vInf*1e3)**2 * params.CD * params.A / 1e3
    
    
    hvec_N = np.cross(xvec_N, vInfvec_N)
    Lupvec_N = np.cross(vInfvec_N, hvec_N)
    
    LupvecU_N = Lupvec_N / np.linalg.norm(Lupvec_N)
    
    vInfvecU_N = vInfvec_N / np.linalg.norm(vInfvec_N)
    
    LvecU_N = LupvecU_N * np.cos(params.bank) + \
              np.cross(vInfvecU_N, LupvecU_N) * np.sin(params.bank)
              
    DvecU_N = - vInfvecU_N
    
    FLvec_N = Lmag * LvecU_N
    FDvec_N = Dmag * DvecU_N

    
    dydt[0] = dx
    dydt[1] = dy
    dydt[2] = dz
    dydt[3:6] = (Fgvec_N + FLvec_N + FDvec_N) / params.m
    
    
    return dydt


def below_min_alt(t, y, params):    
    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    h = r - params.p.rad
    return h - params.hmin

def above_max_alt(t, y, params):
    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    h = r - params.p.rad
    return h - params.hmax