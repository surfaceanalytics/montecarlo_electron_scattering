# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 17:13:41 2020

@author: Mark
"""
import os
import pathlib
os.chdir(pathlib.Path(__file__).parent)

import numpy as np
import matplotlib.pyplot as plt

from montecarlo.simulation import GasPhaseSimulation
from analysis import Analysis

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
#%%
if __name__ == '__main__':
    loss_function = 'He_loss_fn.csv'
    #loss_function=[[10],[0.5]]
    distances = [100000, 300000, 400000, 500000, 800000,900000]
    diameters = [100000, 300000, 400000, 500000, 800000,900000]
    pressures = [0.1, 1,10,25,50]
    nozzles = [100000, 300000, 400000, 500000, 800000,900000]
    angles = [5,10,15,20,30,35,40,50,60]
    Results = []
    counter = 0
    for distance in distances:
        for diameter in diameters:
            for pressure in pressures:
                sim = GasPhaseSimulation(sample_nozzle_distance = distance,
                                      pressure = pressure,
                                      source_diameter = diameter,
                                      initial_KE = 1100,
                                      loss_function = loss_function)
                sim.simulateMany(5000, 'start finish')
                results = sim.start_finish
                A = Analysis(results, parameters = sim.parameters)
                for nozzle in nozzles:
                    nozzle_diameter = nozzle
                    for angle in angles:
                        A.pathHistogram(bins=200, x_limits=(290000,350000),
                                        accept_angle=angle, nozzle_diameter=nozzle_diameter)
                        A.areasUnderProfiles()
                        params = A.parameters
                        p = A.profiles
                        a = A.areas
                        avg = A.getAveragePathLength()
                        params['avg_path'] = avg
                        Results += [params.copy()]
                        print('counter : ' + str(counter))
                        counter+=1

    R_array = []

    keys = ['norm_path','accept_angle','nozzle_diameter',
            'pressure','source_diameter', 'sample_nozzle_distance']
    for R in Results:
        factor = 1000000
        R['norm_path'] = R['avg_path'] / R['sample_nozzle_distance']
        R['nozzle_diameter'] = R['nozzle_diameter'] / factor
        R['source_diameter'] = R['source_diameter'] / factor
        R['sample_nozzle_distance'] = R['sample_nozzle_distance'] / factor
        R_array += [[R[k] for k in keys]]

    R_array = np.array(R_array).copy()
    R_array = R_array[1:,:]

    Y = R_array[:,0]

    X = R_array[:,1:]

    model = Pipeline([('poly', PolynomialFeatures(degree=3)),
                      ('linear', LinearRegression(fit_intercept=False))])

    model = model.fit(X, Y)
    powers = model.steps[0][1].powers_
    coef = model.named_steps['linear'].coef_

    variables = ['x1','x2','x3','x4','x5']

    features = []
    for po in powers:
        f = ''
        for idx, p in enumerate(po):
            if p > 0:
                f += variables[idx] + '^' + str(p)
        features += [f]

    little_x = np.array([[30, 0.1,13, 0.600000,0.200000]])
    avg_predict = model.predict(little_x)[0] * 200000

    trend = []
    for i in range(10):
        z = i / 10
        x = np.array([[22, 1,13, 1,z]])
        avg_predict = model.predict(x)[0]
        trend += [avg_predict]

    plt.plot(trend)