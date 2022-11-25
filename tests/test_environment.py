# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 16:57:44 2021

@author: Mark
"""
import os
import pathlib
os.chdir(pathlib.Path(__file__).parent.parent.parent)

from montecarlo.simulation import GasPhaseSimulation
from montecarlo.analysis import Analysis

#%% Run the simulation (this could take a while).

loss_function = 'He_loss_fn.csv'
sim = GasPhaseSimulation(sample_nozzle_distance = 500000,
                      pressure = 10,
                      source_diameter = 900000,
                      initial_KE = 1100,
                      loss_function = loss_function)
sim.simulateMany(1000, 'start finish')
results = sim.start_finish

print(len(results))

print(sim.scatterer.getIMFP(1100))

#%% Do the analysis

# First create an Analysis object, and pass the results to it
A = Analysis(results, parameters = sim.parameters)
nozzle_diameter = 600000

# Then create a histogram of path lengths. There is one histogram for
# for each n, where n is the number of times inelastically scattered
# The agruments here are: 'show' shows the plot, 'bin' is the number of
# bins in the histogram, x_limits is a tuple that defines the limits of
# the plot, accept_angle takes a tuple, where the first position is angle
# and the second is the radius of the simulation boundary
A.pathHistogram('show', bins=200,
                x_limits=(A.parameters['sample_nozzle_distance'],
                          A.parameters['sample_nozzle_distance']+100000),
                accept_angle=30, nozzle_diameter=nozzle_diameter)

# Then get the areas under the histogram profiles. This is the same as
# counting the number of electrons that have been scatterded n times.
A.areasUnderProfiles('show')

# Then plot the start positions and end positions of each electron
A.plotStartFinish(x_lim = 1000000)

A.showSpectrum(step_size = 0.25, collected_only = False,
               energy_range = (1040,1103))

params = A.parameters
p = A.profiles
a = A.areas
avg = A.getAveragePathLength()

R = A.results
S = A.selection