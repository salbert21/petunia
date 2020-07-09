# -*- coding: utf-8 -*-
"""
Created on Sat May  2 11:53:51 2020

@author: Samuel
"""

import numpy as np
from scipy.integrate import solve_ivp, trapz
import sys

import ODEs
from atm import getRho_from_table
from conversions import LLAEHV2RV, RV2LLAEHV, VN2Vinf, Vinf2VN
import constants

class Params:
       
    class p:
        pass 
    
class Outs:
    pass


def simRun(params, tspan, events, **options):
    # now do the setup and call dynamics
    y0 = np.block([params.x0, params.v0])
        
    sol = solve_ivp(lambda t, y: ODEs.dynamics(t,y,params), 
                    [tspan[0], tspan[-1]], y0.flatten(),
                    rtol=1e-10,atol=1e-10,
                    events=events) 
                    # 1e-10 maintains about 5 sig figs of accuracy
    if 'verbose' in options:
        if options['verbose']:
            print('BC: %.2f, EFPA: %.2f' % (params.BC, params.efpa))
            print(sol.message)
    
    if not sol.success:
        sys.exit('integration failed')
    
    if not sol.status:
        print('Got stuck at:')
        print('BC: %.2f' % params.BC)
        print('EFPA: %.2f' % params.efpa)
        sys.exit('no termination event reached') 
    
    return sol

def mainAD(params, tspan, events, outs):
    ### GET DERIVED VEHICLE PARAMETERS
    # get A from CD and BC
    params.A = params.m / (params.BC * params.CD)
    
    # get Rn from A, assuming Rn/Rb = 1/2
    params.Rn = np.sqrt(params.A / np.pi) / 2
    
    # get CL from L/D and CD
    params.CL = params.CD * params.LD
    
    ### CONVERT GIVEN INPUT TO INERTIAL VECTORS (and other input types as well)
    # NOTE - no matter what the input type is, vectors should be in N frame
    if params.inputType == 'inertial vectors':
        params.vInfvec_N = VN2Vinf(params.x0, params.v0, params, tspan[0])
        params.lat, params.lon, params.alt, params.efpa, params.hda,\
            params.vmag = RV2LLAEHV(params.x0, params.v0, params, tspan[0])
        _, _, _, params.efpaWR, params.hdaWR, params.vmagWR = \
            RV2LLAEHV(params.x0, params.vInfvec_N, params, tspan[0])
    
    elif params.inputType == 'wind-relative vectors':
        params.v0 = Vinf2VN(params.x0, params.vInfvec_N, params, tspan[0])
        params.lat, params.lon, params.alt, params.efpa, params.hda,\
            params.vmag = RV2LLAEHV(params.x0, params.v0, params, tspan[0])
        _, _, _, params.efpaWR, params.hdaWR, params.vmagWR = \
            RV2LLAEHV(params.x0, params.vInfvec_N, params, tspan[0])
        
    elif params.inputType == 'inertial angles':
        params.x0, params.v0 = LLAEHV2RV(params.lat, params.lon, params.alt,
                               params.efpa, params.hda, params.vmag, params,
                               tspan[0])
        params.vInfvec_N = VN2Vinf(params.x0, params.v0, params, tspan[0])
        _, _, _, params.efpaWR, params.hdaWR, params.vmagWR = \
            RV2LLAEHV(params.x0, params.vInfvec_N, params, tspan[0])
        
    elif params.inputType == 'wind-relative angles':
        params.x0, params.vInfvec_N = LLAEHV2RV(params.lat, params.lon,
                                                params.alt, params.efpaWR,
                                                params.hdaWR, params.vmagWR,
                                                params, tspan[0])
        params.v0 = Vinf2VN(params.x0, params.vInfvec_N, params, tspan[0])
        _, _, _, params.efpa, params.hda, params.vmag = \
            RV2LLAEHV(params.x0, params.v0, params, tspan[0])
        
    else:
        sys.exit('input type not recognized')
    
    ### CALL SIMRUN
    sol = simRun(params, tspan, events, verbose=True)
    rvec_N = sol.y[0:3,:]
    vvec_N = sol.y[3:6,:] 
    
    ### GET OUTPUT PARAMS OF INTEREST
    # final inertial state
    rfvec_N = rvec_N[:,-1]
    vfvec_N = vvec_N[:,-1]
    tf = sol.t[-1]
    
    # convert final state params (these are inertial)
    lat, lon, alt, fpa, hda, vmag = RV2LLAEHV(rfvec_N, vfvec_N, params, tf)
    
    ## Compute peak heat rate, add to output
    # get airspeed at each time, vInf
    vInfvec_N = []
    for i in range(vvec_N.shape[1]):
        vInfvec_N.append(VN2Vinf(rvec_N[:,i], vvec_N[:,i], params, sol.t[i]))
    
    vInfvec_N = np.asarray(vInfvec_N).T
    vInf = np.linalg.norm(vInfvec_N, axis=0)
    
    # get density at each time
    r = np.linalg.norm(rvec_N, axis=0)
    h = r - params.p.rad
    if params.dMode == 'fun':
        rho = params.dFun(h)
    elif params.dMode == 'table':
        rho = getRho_from_table(params.atmdat, h)
    else:
        sys.exit('atm mode not recognized')
        
    # calculate S-G heat rate without coefficient, SG
    SG = []
    for i in range(len(rho)):
        SG.append(np.sqrt(rho[i] / params.Rn) * vInf[i]**3)
    
    SG = np.asarray(SG)
    SGpeak = SG.max()
    
    # calculate peak heat rate
    qpeak = params.p.k * SGpeak * 1e5 # puts q in W/cm^2 units
    q = params.p.k * SG * 1e5
    
    # now integrate numerically to get heat load
    Qload = trapz(q, sol.t) # J / cm^2
    
    ## Compute max g, add to output
    # run through dynamics again to get forces output
    Fgvec_N = np.empty((3, len(sol.t)))
    Fgvec_N[:] = np.NaN
    FLvec_N = np.empty((3, len(sol.t)))
    FLvec_N[:] = np.NaN
    FDvec_N = np.empty((3, len(sol.t)))
    FDvec_N[:] = np.NaN
    
    for i, (ti, xi, vi) in enumerate(zip(sol.t, rvec_N.T, vvec_N.T)):
        yyi = np.block([xi, vi])
        Fgvec_N[:,i], FLvec_N[:,i], FDvec_N[:,i] = ODEs.dynamics(ti, yyi,
                                                  params, returnForces=True)
    # compute g-load at each time
    gload = (np.linalg.norm(FLvec_N + FDvec_N - Fgvec_N, axis=0) / params.m)\
        / (constants.G0 / 1e3) # divide g0 by 1000 to get it in km/s^2 units
        
    gpeak = gload.max()
    
    ## Compute apoapsis at final time, add to output
    vf = np.linalg.norm(vfvec_N)
    rf = np.linalg.norm(rfvec_N)
    engf = vf**2 / 2 - params.p.mu / rf
    hfvec_N = np.cross(rfvec_N, vfvec_N)
    hf = np.linalg.norm(hfvec_N)
    
    af = - params.p.mu / (2 * engf)
    eccf = np.sqrt(1 + 2 * engf * hf**2 / params.p.mu**2)
    
    raf = af * (1 + eccf)
    haf = raf - params.p.rad
    
    ### ASSIGN OUTPUTS TO outs CLASS   
    # initial state
    outs.lat0 = params.lat
    outs.lon0 = params.lon
    outs.alt0 = params.alt
    outs.efpa0 = params.efpa
    outs.hda0 = params.hda
    outs.vmag0 = params.vmag
    outs.efpaWR0 = params.efpaWR
    outs.hdaWR0 = params.hdaWR
    outs.vmagWR0 = params.vmagWR
    outs.rvec_N0 = params.x0
    outs.vvec_N0 = params.v0
    outs.vInfvec_N0 = params.vInfvec_N
    
    outs.tspan = tspan
    
    # vehicle
    outs.m = params.m
    outs.A = params.A
    outs.CL = params.CL
    outs.CD = params.CD
    outs.BC = params.BC
    outs.Rn = params.Rn
    outs.bank = params.bank
    
    
    # final state
    outs.latf = lat
    outs.lonf = lon
    outs.altf = alt
    outs.fpaf = fpa
    outs.hdaf = hda
    outs.vmagf = vmag
    outs.rvec_Nf = rfvec_N
    outs.vvec_Nf = vfvec_N
    outs.raf = raf
    outs.haf = haf
    outs.engf = engf
    outs.af = af
    outs.eccf = eccf
    outs.t = tf
    
    # peak values and total loads
    outs.SGpeak = SGpeak
    outs.qpeak = qpeak
    outs.Qload = Qload
    outs.gpeak = gpeak

    
    
    if params.returnTimeVectors:
        outs.rvec_N = rvec_N
        outs.vvec_N = vvec_N
        outs.tvec = sol.t
        outs.q = q
        outs.gload = gload

    return outs