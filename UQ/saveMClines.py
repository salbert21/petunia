# -*- coding: utf-8 -*-
"""
saveMClines.py:
    after running analyzeConvergence, set user input values and run this
Created on Thu Nov 19 21:32:34 2020

@author: Samuel Albert
"""


import numpy as np

data = np.load('convergenceData.npz', allow_pickle = True)
fpafMeanErrListList = data['fpafMeanErrListList']
fpafVarErrListList = data['fpafVarErrListList']
engfMeanErrListList = data['engfMeanErrListList']
engfVarErrListList = data['engfVarErrListList']


# =============================================================================
# Set the average, best, and worst case indicies for each plot. ASSUME p = 2
# =============================================================================

## FPAF MEAN
avg = 2
best = 4
worst = 13

fpafMeanErrAvg = fpafMeanErrListList[:,avg]
fpafMeanErrBest = fpafMeanErrListList[:,best]
fpafMeanErrWorst = fpafMeanErrListList[:,worst]

## FPAF VARIANCE
avg = 12
best = 4
worst = 3

fpafVarErrAvg = fpafVarErrListList[:,avg]
fpafVarErrBest = fpafVarErrListList[:,best]
fpafVarErrWorst = fpafVarErrListList[:,worst]

## ENGF MEAN
avg = 1
best = 7
worst = 13

engfMeanErrAvg = engfMeanErrListList[:,avg]
engfMeanErrBest = engfMeanErrListList[:,best]
engfMeanErrWorst = engfMeanErrListList[:,worst]

## ENGF VARIANCE
avg = 9
best = 8
worst = 3

engfVarErrAvg = engfVarErrListList[:,avg]
engfVarErrBest = engfVarErrListList[:,best]
engfVarErrWorst = engfVarErrListList[:,worst]

# =============================================================================
# Save only these trendlines
# =============================================================================

np.savez('MCtrendlines_p2_70000.npz',
         fpafMeanErrAvg = fpafMeanErrAvg,
         fpafMeanErrBest = fpafMeanErrBest,
         fpafMeanErrWorst = fpafMeanErrWorst,
         fpafVarErrAvg = fpafVarErrAvg,
         fpafVarErrBest = fpafVarErrBest,
         fpafVarErrWorst = fpafVarErrWorst,
         engfMeanErrAvg = engfMeanErrAvg,
         engfMeanErrBest = engfMeanErrBest,
         engfMeanErrWorst = engfMeanErrWorst,
         engfVarErrAvg = engfVarErrAvg,
         engfVarErrBest = engfVarErrBest,
         engfVarErrWorst = engfVarErrWorst)


















