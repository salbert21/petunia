# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:52:56 2020

@author: Samuel
"""

import numpy as np
from frames import DCMs
import sys

from atm import getRho_from_table, getWind

def cross(v1,v2):
    '''
    cross product function to speed up computation. Must take exactly 2 3-vecs
    '''
    s1 = v1[1] * v2[2] - v1[2] * v2[1]
    s2 = v1[2] * v2[0] - v1[0] * v2[2]
    s3 = v1[0] * v2[1] - v1[1] * v2[0]
    return np.array([s1, s2, s3])

def norm(v):
    '''
    vector norm function to speed up computation. Must take exactly 1 3-vec
    '''
    return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def sphericalEntryEOMsAug(t, yy, sig0, e0, params):
    '''
    Like sphericalEntryEOMs from ASEN 6015, but dimensional and with GRAM dens.
    Defines 3D 3DOF EOMs for lifting entry vehicle over *ellipsoidal* rotating
        planet. Defined as in FNPAG (Lu et. al), then augmented with 7th state
        to capture evolution of range. Uses bank angle profile as defined in
        FNPEG.
        ASSUMPTIONS:
            - constant bank angle
            - exponential atmosphere
            - rotating planet
            - ellipsoidal planet (J2)
    params requirements:
        OmP: rotation rate of planet, rad/s
        radP: equatorial radius of planet, m
        mu: gravitational parameter of planet, m^3/s^2
        J2: zonal coefficient J2 of planet
        L_D: vehicle L/D ratio
        BC: vehicle ballistic coefficient, kg/m^2
    INPUTS:
        t: time (not used)
        yy: state vector:
            r: radial distance, m
            lon: longitude, radians
            lat: latitude, radians
            v: planet-relative velocity magnitude, m/s
            gam: flight path angle, negative-down, radians
            hda: heading angle, clockwise from north, radians
        sig: bank angle, rad
    OUTPUTS:
        dydt: time derivative of each state
    '''
    # extract constants from params
    OmP = params.p.om
    radP = params.p.rad * 1e3
    mu = params.p.mu * 1e9 # convert from km^3/s^2 to m^3/s^2
    J2 = params.p.J2
    L_D = params.LD
    BC = params.BC
    
    ef = params.ef
    sigf = params.sigf
    
    # extract state variables
    r = yy[0]
    # lon = yy[1] # not used in EOMs
    lat = yy[2]
    v = yy[3]
    gam = yy[4]
    hda = yy[5]
    
    # get gravity terms from J2 model
    gr = mu/r**2 * (1 + J2 * (radP/r)**2 * (1.5 - 4.5 * np.sin(lat)**2))
    gphi = mu/r**2 * (J2 * (radP/r)**2 * (3 * np.sin(lat) * np.cos(lat)))
    
    # get density at this altitude
    if params.dMode == 'fun':
        rho = params.dFun((r - radP)/1e3)
    elif params.dMode == 'table':
        rho = getRho_from_table(params.atmdat, (r - radP)/1e3)
    else:
        sys.exit('atm mode not recognized')
        
    # #### TROUBLESHOOTING CODE ####
    # rho0 = 0.0263 # kg/m^3
    # H = 10153 # m
    # h = r - radP
    # rho = rho0 * np.exp(-h / H)
    # #### TROUBLESHOOTING CODE ####
    
    # compute e
    e = mu / r - v**2 / 2
    
    # compute current bank angle
    sig = sig0 + (e - e0) / (ef - e0) * (sigf - sig0)
    sig = np.radians(sig)
    
    # compute lift and drag accelerations
    D = rho * v**2 / (2 * BC)
    L = L_D * D
    
    # EOMs
    dydt = np.empty(7)
    dydt[:] = np.NaN
    
    dydt[0] = v * np.sin(gam)
    dydt[1] = v * np.cos(gam) * np.sin(hda) / (r * np.cos(lat))
    dydt[2] = v * np.cos(gam) * np.cos(hda) / r
    
    dydt[3] = -D\
              - gr * np.sin(gam)\
              - gphi * np.cos(gam) * np.cos(hda)\
              + OmP**2 * r * np.cos(lat) * (np.sin(gam) * np.cos(lat)\
                                            - np.cos(gam) * np.sin(lat)\
                                                * np.cos(hda))
    dydt[4] = 1/v * (L * np.cos(sig)\
                     + (v**2/r - gr) * np.cos(gam)\
                     + gphi * np.sin(gam) * np.cos(hda)\
                     + 2 * OmP * v * np.cos(lat) * np.sin(hda)\
                     + OmP**2 * r * np.cos(lat) * (np.cos(gam) * np.cos(lat)\
                                                   + np.sin(gam) * np.cos(hda)\
                                                       * np.sin(lat)))
    dydt[5] = 1/v * (L * np.sin(sig) / np.cos(gam)\
                     + v**2/r * np.cos(gam) * np.sin(hda) * np.tan(lat)\
                     + gphi * np.sin(hda) / np.cos(gam)\
                     - 2 * OmP * v * (np.tan(gam) * np.cos(hda) * np.cos(lat)\
                                      - np.sin(lat))\
                     + OmP**2 * r / np.cos(gam) * np.sin(hda)\
                         * np.sin(lat) * np.cos(lat))
    dydt[6] =  - v * np.cos(gam) / r
    
    return dydt

def sphericalEntryEOMs(t, yy, sig, params):
    '''
    defines 3D 3DOF EOMs for lifting entry vehicle over *ellipsoidal* rotating
        planet. Defined as in FNPAG (Lu et. al).
        ASSUMPTIONS:
            - constant bank angle
            - exponential atmosphere
            - rotating planet
            - ellipsoidal planet (J2)
    params requirements:
        OmP: rotation rate of planet, rad/s
        radP: equatorial radius of planet, m
        mu: gravitational parameter of planet, m^3/s^2
        J2: zonal coefficient J2 of planet
        L_D: vehicle L/D ratio
        BC: vehicle ballistic coefficient, kg/m^2
    INPUTS:
        t: time (not used)
        yy: state vector:
            r: radial distance, m
            lon: longitude, radians
            lat: latitude, radians
            v: planet-relative velocity magnitude, m/s
            gam: flight path angle, negative-down, radians
            hda: heading angle, clockwise from north, radians
        sig: bank angle, rad
    OUTPUTS:
        dydt: time derivative of each state
    '''
    # extract constants from params
    OmP = params.p.om
    radP = params.p.rad * 1e3
    mu = params.p.mu * 1e9 # convert from km^3/s^2 to m^3/s^2
    J2 = params.p.J2
    L_D = params.LD
    BC = params.BC
    
    # extract state variables
    r = yy[0]
    # lon = yy[1] # not used in EOMs
    lat = yy[2]
    v = yy[3]
    gam = yy[4]
    hda = yy[5]
    
    # get gravity terms from J2 model
    gr = mu/r**2 * (1 + J2 * (radP/r)**2 * (1.5 - 4.5 * np.sin(lat)**2))
    gphi = mu/r**2 * (J2 * (radP/r)**2 * (3 * np.sin(lat) * np.cos(lat)))
    
    # get density at this altitude
    if params.dMode == 'fun':
        rho = params.dFun((r - radP)/1e3)
    elif params.dMode == 'table':
        rho = getRho_from_table(params.atmdat, (r - radP)/1e3)
    else:
        sys.exit('atm mode not recognized')
        
    # #### TROUBLESHOOTING CODE ####
    # rho0 = 0.0263 # kg/m^3
    # H = 10153 # m
    # h = r - radP
    # rho = rho0 * np.exp(-h / H)
    # #### TROUBLESHOOTING CODE ####
    
    # compute lift and drag accelerations
    D = rho * v**2 / (2 * BC)
    L = L_D * D
    
    # EOMs
    dydt = np.empty(6)
    dydt[:] = np.NaN
    
    dydt[0] = v * np.sin(gam)
    dydt[1] = v * np.cos(gam) * np.sin(hda) / (r * np.cos(lat))
    dydt[2] = v * np.cos(gam) * np.cos(hda) / r
    
    dydt[3] = -D\
              - gr * np.sin(gam)\
              - gphi * np.cos(gam) * np.cos(hda)\
              + OmP**2 * r * np.cos(lat) * (np.sin(gam) * np.cos(lat)\
                                            - np.cos(gam) * np.sin(lat)\
                                                * np.cos(hda))
    dydt[4] = 1/v * (L * np.cos(sig)\
                     + (v**2/r - gr) * np.cos(gam)\
                     + gphi * np.sin(gam) * np.cos(hda)\
                     + 2 * OmP * v * np.cos(lat) * np.sin(hda)\
                     + OmP**2 * r * np.cos(lat) * (np.cos(gam) * np.cos(lat)\
                                                   + np.sin(gam) * np.cos(hda)\
                                                       * np.sin(lat)))
    dydt[5] = 1/v * (L * np.sin(sig) / np.cos(gam)\
                     + v**2/r * np.cos(gam) * np.sin(hda) * np.tan(lat)\
                     + gphi * np.sin(hda) / np.cos(gam)\
                     - 2 * OmP * v * (np.tan(gam) * np.cos(hda) * np.cos(lat)\
                                      - np.sin(lat))\
                     + OmP**2 * r / np.cos(gam) * np.sin(hda)\
                         * np.sin(lat) * np.cos(lat))
    
    return dydt

def dynamics(t, yy, sig, params, **options):
    '''
    ODEs for full dynamics acting on the vehicle
    assume:
        - just one vehicle
        - yy is [position; velocity] in INERTIAL frame
        - sig is bank angle in deg
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
    
    r = norm(xvec_N)
    mu_r3 = params.p.mu / r**3
    
    ## Get gravitational force
    Fgvec_N = - mu_r3 * params.m * xvec_N
    
    ## Get aerodynamics forces
    # TODO: carefully replace below with VN2Vinf function call
    SN = DCMs.getSN(t, params)
    NS = DCMs.getNS(t, params)
    
    OMvec = np.array([0, 0, -params.p.om])
    
    vvec_S = SN @ (vvec_N + cross(OMvec, xvec_N))
    
    h = r - params.p.rad
    
    if params.dMode == 'fun':
        rho = params.dFun(h)
    elif params.dMode == 'table':
        rho = getRho_from_table(params.atmdat, h)
    else:
        sys.exit('atm mode not recognized')
    wvec_S = getWind(h)
    
    vInfvec_S = vvec_S + wvec_S
    vInf = norm(vInfvec_S)
    
    vInfvec_N = NS @ vInfvec_S
    
    Lmag = 1/2 * rho * (vInf*1e3)**2 * params.CL * params.A / 1e3
    Dmag = 1/2 * rho * (vInf*1e3)**2 * params.CD * params.A / 1e3
    
    
    hvec_N = cross(xvec_N, vInfvec_N)
    Lupvec_N = cross(vInfvec_N, hvec_N)
    
    LupvecU_N = Lupvec_N / norm(Lupvec_N)
    
    # TODO - following code check is a little sketchy. can improve w/ 6 DOF
    c1 = cross(xvec_N, params.v0)
    c2 = cross(xvec_N, vInfvec_N)
    # if we are "flying upside down"
    if np.dot(c1,c2) < 0:
        # then change the sign on the lift vector to be upside down
        LupvecU_N = -LupvecU_N
    
    # print(LupvecU_N)
    # print(t)
    
    vInfvecU_N = vInfvec_N / norm(vInfvec_N)
    
    LvecU_N = LupvecU_N * np.cos(np.radians(sig)) + \
              cross(vInfvecU_N, LupvecU_N) * np.sin(np.radians(sig))
              
    DvecU_N = - vInfvecU_N
    
    FLvec_N = Lmag * LvecU_N
    FDvec_N = Dmag * DvecU_N

    
    dydt[0] = dx
    dydt[1] = dy
    dydt[2] = dz
    dydt[3:6] = (Fgvec_N + FLvec_N + FDvec_N) / params.m
    
    if 'returnForces' in options:
        if options['returnForces']:
            return Fgvec_N, FLvec_N, FDvec_N
    else:
        return dydt


def below_min_alt(t, y, params):    
    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    h = r - params.p.rad
    return h - params.hmin

def above_max_alt(t, y, params):
    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    h = r - params.p.rad
    return h - params.hmax

def below_min_alt_sph(t, y, params):    
    r = y[0]/1e3
    h = r - params.p.rad
    return h - params.hmin

def above_max_alt_sph(t, y, params):
    r = y[0]/1e3
    h = r - params.p.rad
    return h - params.hmax

def switchEvent(t, y, ts):
    return t - ts