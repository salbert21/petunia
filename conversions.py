# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 10:10:41 2020

@author: Samuel Albert
"""

import numpy as np
import scipy.linalg as LA
import sys

from atm import getWind
from frames import DCMs

def LLA2R(lat, lon, alt, rad):
    r = rad + alt
    x = r * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
    y = r * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
    z = r * np.sin(np.radians(lat))
    
    rvec_S = np.array([x,y,z])
    return rvec_S

def LLAEHV2RV(lat, lon, alt, fpa, hda, vmag, params, t):
    '''
    takes lat, long, alt, fpa, heading angle, and velocity magnitude
    returns position and velocity vectors in inertial frame
    '''
    
    # get position vector in surface frame components
    rvec_S = LLA2R(lat, lon, alt, params.p.rad)
    
    # get NED unit vectors in surface frame components
    ## first assume lon = 0:
    uD_S = -np.cos(np.radians(lat)) * np.array([1,0,0]) + \
           -np.sin(np.radians(lat)) * np.array([0,0,1])
    uN_S = -np.sin(np.radians(lat)) * np.array([1,0,0]) + \
            np.cos(np.radians(lat)) * np.array([0,0,1])
    uE_S = np.array([0,1,0])
    
    ## now rotate by longitude:
    M3T = DCMs.getM3(np.radians(lon)).T
    uN_S = M3T @ uN_S
    uE_S = M3T @ uE_S
    uD_S = M3T @ uD_S
    
    # find vvec_S using NED vectors
    vD_S = vmag * -np.sin(np.radians(fpa))
    vH_S = vmag * np.cos(np.radians(fpa))
    vE_S = vH_S * np.cos(np.radians(hda))
    vN_S = vH_S * -np.sin(np.radians(hda))
    
    vvec_S = vN_S * uN_S + vE_S * uE_S + vD_S * uD_S
    
    # rotate vectors to inertial frame
    #### NOTE: this is only a rotation, not a transport ####
    NS = DCMs.getNS(t, params)
    rvec_N = NS @ rvec_S
    vvec_N = NS @ vvec_S
    
    return rvec_N, vvec_N

def RV2LLAEHV(rvec_N, vvec_N, params, t):
    
    # rotate vectors to surface frame
    SN = DCMs.getSN(t, params)
    rvec_S = SN @ rvec_N
    vvec_S = SN @ vvec_N
    
    # get lat, lon, alt
    r = LA.norm(rvec_S)
    alt = r - params.p.rad
    
    lat = np.degrees(np.arcsin(rvec_S[2]/r))
    lon = np.degrees(np.arctan2(rvec_S[1], rvec_S[0]))
    
    # get vmag and fpa
    vmag = LA.norm(vvec_S)
    fpa = np.degrees(np.arcsin(np.dot(rvec_S,vvec_S)\
                               / (LA.norm(rvec_S)*LA.norm(vvec_S))))
    
    # now get hda
    # get NED unit vectors in surface frame components
    ## first assume lon = 0:
    uD_S = -np.cos(np.radians(lat)) * np.array([1,0,0]) + \
           -np.sin(np.radians(lat)) * np.array([0,0,1])
    uN_S = -np.sin(np.radians(lat)) * np.array([1,0,0]) + \
            np.cos(np.radians(lat)) * np.array([0,0,1])
    uE_S = np.array([0,1,0])
    
    ## now rotate by longitude:
    M3T = DCMs.getM3(np.radians(lon)).T
    uN_S = M3T @ uN_S
    uE_S = M3T @ uE_S
    uD_S = M3T @ uD_S
    
    ## now project vvec_S onto horizontal plane by subtracting uD_S component
    vvecH_S = vvec_S - np.dot(vvec_S,uD_S) / (LA.norm(uD_S)**2) * uD_S
    
    ## finally, find angle between the horizontal projection of vvec and East
    # check to make sure dot product doesn't out of range numerically
    val = np.dot(vvecH_S,uE_S) / (LA.norm(vvecH_S) * LA.norm(uE_S))
    ep = 1e-5
    if val > 1:
        if val-ep > 1:
            sys.exit('value error, arccos of value larger than 1')
        else:
            val = 1
    elif val < -1:
        if val+ep < -1:
            sys.exit('value error, arccos of value less than -1')
        else:
            val = -1
            
    hda = np.degrees(np.arccos(val))
        
    # make sure hda has the right sign
    if np.dot(uN_S, vvecH_S) > 0: # if vvec points up toward north
        hda = -hda # flip the sign on the angle
        
    return lat, lon, alt, fpa, hda, vmag
    
    
def VN2Vinf(rvec_N, vvec_N, params, t):
    '''
    converts an inertial velocity vector to a wind-relative velocity vector
    (airspeed vector), returned in INTERTIAL FRAME COMPONENTS
    '''
    
    # get DCMs
    SN = DCMs.getSN(t, params)
    NS = DCMs.getNS(t, params)
    
    # get rotation vector of inertial frame w.r.t. surface frame, om N/S, 
    #  hence the negative sign on om
    OMvec = np.array([0, 0, -params.p.om])
    
    # transport and rotate to get groundspeed vec in surface components
    vvec_S = SN @ (vvec_N + np.cross(OMvec, rvec_N))
    
    # get wind in surface frame components (dummy function currently)
    r = np.linalg.norm(rvec_N)
    h = r - params.p.rad
    wvec_S = getWind(h)
    
    # sum to get total airspeed vector in surface frame
    vInfvec_S = vvec_S + wvec_S
    
    # now rotate airspeed into inertial frame components and return
    vInfvec_N = NS @ vInfvec_S
    
    return vInfvec_N

def Vinf2VN(rvec_N, vInfvec_N, params, t):
    '''
    converts an airspeed vector in inertial frame to an inertial velocity
    vector in inertial components
    '''
    
    # get DCMs
    SN = DCMs.getSN(t, params)
    NS = DCMs.getNS(t, params)
    
    rvec_S = SN @ rvec_N
    
    # rotate airspeed vec to surface frame
    vInfvec_S = SN @ vInfvec_N
    
    # get wind and subtract from airspeed to get groundspeed
    r = np.linalg.norm(rvec_N)
    h = r - params.p.rad
    wvec_S = getWind(h)
    
    vvec_S = vInfvec_S - wvec_S
    
    # transport and rotate to get inertial velocity in inertial frame:
        
    # get rotation vector of surface w.r.t. inertial frame, om_S/N, 
    OMvec = np.array([0, 0, params.p.om])
    
    vvec_N = NS @ (vvec_S + np.cross(OMvec, rvec_S))
    
    return vvec_N


def getApses(rvec_N, vvec_N, params):
    '''
    converts inertial cartesian vectors to orbital energy and apoapsis altitude
    '''
    
    r = np.linalg.norm(rvec_N)
    v = np.linalg.norm(vvec_N)
    
    eng = v**2/2 - params.p.mu / r
    
    hvec_N = np.cross(rvec_N, vvec_N)
    h = np.linalg.norm(hvec_N)
    
    a = - params.p.mu / (2 * eng)
    e = np.sqrt(1 + 2 * eng * h**2 / params.p.mu**2)
    
    ra = a * (1 + e)
    rp = a * (1 - e)
    
    return ra, rp

def getApsesSphPR(xxvec, params, returnEng = False):
    '''
    returns radius of periapsis and apoapsis for given spherical state.
    ASSUMES given state is PLANET-RELATIVE
    INPUTS:
        xxvec: planet-relative spherical state, m, m/s, rad
        mu: planet gravitational parameter, m^3/s^2
    OUTPUTS:
        ra: radius of apoapsis, m
        rp: radius of periapsis, m
        ENG: specific mechanical orbital energy, m^2/s^2
    '''
    
    mu = params.p.mu * 1e9
    
    # get inertial cartesian state
    rvec, vrvec = sph2cart(xxvec)
    vivec = vr2viCart(rvec, vrvec, params)
    r = np.linalg.norm(rvec)
    v = np.linalg.norm(vivec)
    
    # get Keplerian parameters
    hvec = np.cross(rvec, vivec)
    h = np.linalg.norm(hvec)
    ENG = v**2/2 - mu/r
    SMA = - mu / (2 * ENG)
    arg = 1 + 2 * ENG * h**2 / mu**2
    if abs(arg) < 1e-12:
        ECC = 0 # circular orbit case
    else:
        ECC = np.sqrt(arg)
    
    # get apses radii
    ra = SMA * (1 + ECC)
    rp = SMA * (1 - ECC)
    
    if returnEng:
        return ra, rp, ENG
    else:
        return ra, rp

def sph2cart(xsphvec):
    '''
    converts spherical to cartesian coordinates for EDL.
    Uses ASEN 6015 derivations.
    INPUTS:
        xphvec:
            r: radius
            lon: longitude, rad
            lat: latitude, rad
            v: velocity magnitude
            gam: flight path angle, negative-down, rad
            hda: heading angle, 0-North clockwise, rad
    OUTPUTS:
        rvec: position vector
        vvec: velocity vector
    '''
    
    # extract state variables
    r = xsphvec[0]
    lon = xsphvec[1]
    lat = xsphvec[2]
    v = xsphvec[3]
    gam = xsphvec[4]
    hda = xsphvec[5]
    
    # get radius vector
    rX = r * np.cos(lat) * np.cos(lon)
    rY = r * np.cos(lat) * np.sin(lon)
    rZ = r * np.sin(lat)
    rvec = np.array([rX, rY, rZ])
    
    # get velocity vector
    vE = v * np.cos(gam) * np.sin(hda)
    vN = v * np.cos(gam) * np.cos(hda)
    vR = v * np.sin(gam)
    
    vX = -vE * np.sin(lon)\
          - vN * np.sin(lat) * np.cos(lon)\
          + vR * np.cos(lat) * np.cos(lon)
    vY = vE * np.cos(lon)\
          - vN * np.sin(lat) * np.sin(lon)\
          + vR * np.cos(lat) * np.sin(lon)
    vZ = vN * np.cos(lat) + vR * np.sin(lat)
    
    vvec = np.array([vX, vY, vZ])
    
    return rvec, vvec

def cart2sph(rvec, vvec):
    '''
    converts cartesian to spherical coordinates for EDL.
    Uses ASEN 6015 derivations.
    INPUTS:
        rvec: position vector
        vvec: velocity vector
    OUTPUTS:
        xsphvec:
            r: radius
            lon: longitude, rad
            lat: latitude, rad
            v: velocity magnitude
            gam: flight path angle, negative-down, rad
            hda: heading angle, 0-North clockwise, rad
    '''
    
    r = np.linalg.norm(rvec)
    v = np.linalg.norm(vvec)
    
    lon = np.arctan2(rvec[1], rvec[0])
    lat = np.arcsin(rvec[2]/r)
    
    vE = -vvec[0] * np.sin(lon) + vvec[1] * np.cos(lon)
    vN = -vvec[0] * np.sin(lat) * np.cos(lon)\
        - vvec[1] * np.sin(lat) * np.sin(lon) + vvec[2] * np.cos(lat)
    hda = np.arctan2(vE, vN)
    
    vR = vvec[0] * np.cos(lat) * np.cos(lon)\
        + vvec[1] * np.cos(lat) * np.sin(lon) + vvec[2] * np.sin(lat)
    gam = np.arcsin(vR / v)
    
    return np.array([r, lon, lat, v, gam, hda])
    

def vr2viCart(rvec, vrvec, params):
    '''
    converts planet-relative velocity to inertial velocity assuming no wind.
    INPUTS:
        rvec: position vector, m
        vvec: planet-relative velocity vector, m/s
        OmP: planet self-rotation rate (in positive z direction), rad/s
    OUTPUTS:
        vivec: inertial velocity vector, m/s
    '''
    
    Omvec = np.array([0, 0, params.p.om])
    vivec = vrvec + np.cross(Omvec, rvec)
    return vivec
    
    
def getAzimuth(lon1, lat1, lon2, lat2):
    '''
    returns bearing (Azimuth angle along great circle) between two points.
    Input angles must be in radians. output is also in radians
    '''
    
    x = np.cos(lat2) * np.sin(lon2 - lon1)
    y = np.cos(lat1) * np.sin(lat2)\
        - np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1)
    Az = np.arctan2(x,y)
    
    return Az

def greatCircleDist(lon1, lat1, lon2, lat2, params):
    '''
    Computes the great-circle distance between two points in km
    '''
    
    dlon = abs(lon2 - lon1)
    dsig = np.arccos(np.sin(lat1) * np.sin(lat2)\
                     + np.cos(lat1) * np.cos(lat2) * np.cos(dlon))
    d = dsig * params.p.rad
    
    return d
    































    
    
    