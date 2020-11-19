# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 21:32:41 2020

@author: Samuel Albert
"""

'''NOTE: assumes all combined files have the same dispersed inputs!'''

import numpy as np

MCfilenameList = [
    './results/' + 'Mars_10000_1118131943.npz',
    './results/' + 'Mars_10000_1118132008.npz',
    './results/' + 'Mars_10000_1118132011.npz',
    './results/' + 'Mars_10000_1118132408.npz',
    './results/' + 'Mars_10000_1118132411.npz',
    './results/' + 'Mars_10000_1118132413.npz',
    './results/' + 'Mars_10000_1118142637.npz'
    ]

paramsList = []
outsList = []
m_YList = []
CD_YList = []
vmag_YList = []
vmag_YList = []
efpa_YList = []
atm_YsList = []

for MCfilename in MCfilenameList:
    data = np.load(MCfilename, allow_pickle = True)
    params = data['paramsList']
    paramsList.extend(params)
    outs = data['outsList']
    outsList.extend(outs)
    m_Y = data['m_YList']
    m_YList.extend(m_Y)
    CD_Y = data['CD_YList']
    CD_YList.extend(CD_Y)
    vmag_Y = data['vmag_YList']
    vmag_YList.extend(vmag_Y)
    efpa_Y = data['efpa_YList']
    efpa_YList.extend(efpa_Y)
    atm_Ys = data['atm_YsList']
    atm_YsList.extend(atm_Ys)
    
    CDmean = data['CDmean']
    CDstd = data['CDstd']
    mmean = data['mmean']
    mstd = data['mstd']
    efpamean = data['efpamean']
    efpastd = data['efpastd']
    vmagmean = data['vmagmean']
    vmagstd = data['vmagstd']
    
outname = './results/' + 'Mars_70000_1118142637.npz'

np.savez(outname,
          paramsList = paramsList,
          outsList = outsList,
          m_YList = m_YList,
          CD_YList = CD_YList,
          vmag_YList = vmag_YList,
          efpa_YList = efpa_YList,
          atm_YsList = atm_YsList,
          CDmean = CDmean,
          CDstd = CDstd,
          mmean = mmean,
          mstd = mstd,
          efpamean = efpamean,
          efpastd = efpastd,
          vmagmean = vmagmean,
          vmagstd = vmagstd
          )
    

