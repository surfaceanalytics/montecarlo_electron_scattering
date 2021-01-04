# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:50:34 2021

@author: Mark
"""

import time

from solid_state_sim import SolidStateSimulation
from analysis import Analysis

#%% Run the simulation (this could take a while).

if __name__ == '__main__':
    loss_function = 'Ag_loss_fn.csv'

    sim = SolidStateSimulation(sample_nozzle_distance = 300000,
                          density = 10.4,
                          molar_mass = 47,
                          source_diameter = 300000,
                          initial_KE = 1100,
                          loss_function = loss_function,
                          source_thickness = 10)
    start = time.time()
    sim.simulateMany(2000, 'start finish')
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
    A.writeExcel('nozz300um, dist300um, Ag, ang90deg, 1Me')
    
