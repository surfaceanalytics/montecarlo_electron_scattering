# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:50:34 2021

@author: pielsticker
"""
import os
import pathlib
os.chdir(pathlib.Path(__file__).parent)

import time

from montecarlo.simulation import SolidStateSimulation
from montecarlo.analysis import Analysis

#%% Run the simulation (this could take a while).

if __name__ == '__main__':
    loss_function = 'Ag_loss_fn.csv'

    sample_nozzle_distance = 500000
    pressure = 10
    source_diameter= 300000
    initial_KE = 1100
    number_of_electrons = 100000

    sim = SolidStateSimulation(
        sample_nozzle_distance=sample_nozzle_distance,
        density = 10.4,
        molar_mass = 47,
        source_diameter=source_diameter,
        initial_KE = 1100,
        loss_function = loss_function,
        source_thickness = 10)
    start = time.time()
    sim.simulateMany(number_of_electrons, 'start finish')
    stop = time.time()


    print(stop-start)
    results = sim.start_finish

#%% Do the analysis
if __name__ == '__main__':
    # First create an Analysis object, and pass the results to it
    A = Analysis(results, parameters = sim.parameters)
    nozzle_diameter = 300000

    # Then create a histogram of path lengths. There is one histogram for
    # for each n, where n is the number of times inelastically scattered
    # The agruments here are: 'show' shows the plot, 'bin' is the number of
    # bins in the histogram, x_limits is a tuple that defines the limits of
    # the plot, accept_angle takes a tuple, where the first position is angle
    # and the second is the radius of the simulation boundary
    A.pathHistogram('show',bins=200, x_limits=(0,50),accept_angle=90, nozzle_diameter=nozzle_diameter)

    # Then get the areas under the histogram profiles. This is the same as
    # counting the number of electrons that have been scatterded n times.
    A.areasUnderProfiles('show')

    # Then plot the start positions and end positions of each electron
    #A.plotStartFinish(x_lim = 1000000)

    A.showSpectrum(step_size = 0.25, collected_only = True, energy_range = (1020,1102))

    params = A.parameters
    p = A.profiles
    a = A.areas
    avg = A.getAveragePathLength()

    R = A.results
    S = A.selection

#%%
    filename = (
        f'nozz{int(nozzle_diameter/1000)}um, '
        f'dist{int(sample_nozzle_distance/1000)}um, '
        f'{loss_function.split("_")[0]}, '
        f'ang{A.parameters["accept_angle"]}deg, '
        f'{number_of_electrons}e'
        )
    A.writeExcel(filename)