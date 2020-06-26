# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 10:10:41 2020

@author: Samuel Albert
"""

import numpy as np
import scipy.linalg as LA
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
    hda = np.degrees(np.arccos(np.dot(vvecH_S,uE_S)\
                                    / (LA.norm(vvecH_S) * LA.norm(uE_S))))
        
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
    
    
    
    
    
    