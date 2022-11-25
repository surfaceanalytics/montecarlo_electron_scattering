# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 10:10:53 2020

@author: Mark
"""
import collections
import multiprocessing
from solid_state_sim import SolidStateSimulation

loss_function = 'Ag_loss_fn.csv'
def runOnce():
    sim = SolidStateSimulation(sample_nozzle_distance = 300000,
                  density = 10.4,
                  molar_mass = 47,
                  source_diameter = 300000,
                  initial_KE = 1100,
                  loss_function = loss_function,
                  source_thickness = 30)
    sim.simulateMany(1000, 'start finish')
    return sim.start_finish


if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=2)

    result = [pool.apply(runOnce) for i in range(2)]
