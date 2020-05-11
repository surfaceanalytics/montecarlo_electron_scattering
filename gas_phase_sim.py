# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:49:49 2020

@author: Mark
"""
from simulation import Simulation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from generate_angle import AngleDist
from matplotlib import cm
import pickle
from shapes import Sphere
from analysis import Analysis
#%%

if __name__ == '__main__':
    sim_radius = 800000
    sim = Simulation(sim_radius) 
    d = sim.densityFromP(10) # this gets the density in atoms/nm^3 from pressure
    # argument is the radius of the spherical simulation evnrionment
    sim.initial_KE = 1100
    sim.addSource(400000,1) # source is a shape that emits electrons
    # arguments are radius and thickness (in nm) of a disc
    sim.addScatteringMedium(Sphere(sim_radius))
    sim.addScatterer(d,0.4,0.85,0.95,1.12,4) 
    # arguments are denstiy in atoms/nm^3
    # inel_factor, and inel_exp which determines the inelastic cross section
    # from a fit to TPP2M and NIST data using sigma = facter / (KE * exp)
    # el_factor, which determines the elastic cross section
    # and atomic number
    
    sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford', 
                                                    energy = sim.initial_KE,
                                                    Z = sim.scatterer.Z,
                                                    param = 0.9), 
                           'inelastic':AngleDist(kind = 'Constant')}
    
    
    '''
    sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Constant'), 
                           'inelastic':AngleDist(kind = 'Constant')}'''
    #set the high of the source slice
    sim.height = 1

#%% 

# potentially very long calculation
if __name__ == '__main__':

    sim.simulateMany(200000, 'start finish')
    results = sim.start_finish

#%%
if __name__ == '__main__':
    # First create an Analysis object, and pass the results to it
    A = Analysis(results)
    # Then create a histogram of path lengths. There is one histogram for 
    # for each n, where n is the number of times inelastically scattered
    # The agruments here are: 'show' shows the plot, 'bin' is the number of
    # bins in the histogram, x_limits is a tuple that defines the limits of
    # the plot, accept_angle takes a tuple, where the first position is angle
    # and the second is the radius of the simulation boundary
    A.pathHistogram('show',bins=25, x_limits=(600000,1100000), accept_angle =
                    (30,800000))

    # Then get the areas under the histogram profiles. This is the same as 
    # counting the number of electrons that have been scatterded n times.        
    A.areasUnderProfiles('show')
    
    # Then plot the start positions and end positions of each electron
    A.plotStartFinish(accept_angle = (30,800000))

    p = A.profiles
    a = A.areas
