# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 10:48:10 2020

@author: Mark
"""
from analysis import Analysis

import pickle
import numpy as np

filenames = ['sim_output0','sim_output1','sim_output2','sim_output3',
             'sim_output4','sim_output5','sim_output6','sim_output7',
             'sim_output8','sim_output9','sim_output10','sim_output11',
             'sim_output12', 'sim_output13','sim_output14','sim_output15', 
             'sim_output16','sim_output17','sim_output18','sim_output19',
             'sim_output20','sim_output21','sim_output22','sim_output23',
             'sim_output24','sim_output25','sim_output26','sim_output27',
             'sim_output28','sim_output29','sim_output30','sim_output31']

results = []
n = 0
for filename in filenames:
    with open(filename, 'rb') as f:
        load = pickle.load(f)
        results += [load['result']]
        params = load['parameters']
        n += params['n']
params['n'] = n    

def saveResultsToPickle(result, params, filename):
    to_pickle = {'result':result.copy(), 'parameters':params.copy()}
    with open(filename, 'wb') as f: 
        pickle.dump(to_pickle, f)
       
def loadResultsFromPickle(filename):
    r = pickle.load(open(filename, "rb" ))
    results = r['result']
    params = r['parameters']
    return results, params        
    
def appendResults(results):
    for jdx, j in enumerate(results):
        if jdx == 0:
            R = j
        else:
            for idx, i in enumerate(R):
                R[idx] = np.append(R[idx],j[idx], axis=0)
    return R

results = appendResults(results)
#results, params = loadResultsFromPickle('25MbarHe')

A = Analysis(results, parameters = params)
nozzle_diameter = 300000    
A.pathHistogram('show',bins=100, x_limits=(0,50),accept_angle=90, nozzle_diameter=nozzle_diameter)    
A.areasUnderProfiles('show')
A.showSpectrum(step_size = 0.25, collected_only = True, energy_range = (0,1110))
A.showPartialSpectra(step_size = 0.25, collected_only = True, energy_range = (1000,1110))

A.writeExcel('nozz300um, dist300um, He, angle_spread, 2Me')
A.writeExcel('nozz300um, dist300um, Ag, all, spherical, 3.2Me')

        
saveResultsToPickle(results, params, 'Ag_spherical_angle')
