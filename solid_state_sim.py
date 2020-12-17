# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 14:38:02 2020

@author: Mark
"""

from simulation import Simulation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from generate_angle import AngleDist
from matplotlib import cm
import pickle
from shapes import Sphere, Disc
from analysis import Analysis
import time
#%%

class solidStateSimulation(Simulation):
    def __init__(self, **kwargs):
        
        '''This will be used as the radius of the simulation boundary'''
        if 'sim_radius' in kwargs.keys():
            sim_radius = kwargs['sim_radius']
        else:
            sim_radius = 1000000
            
        if 'sample_nozzle_distance' in kwargs.keys():
            sample_nozzle_distance = kwargs['sample_nozzle_distance']
        else:
            sample_nozzle_distance = 300000
        
        ''' This is the initial kinetic energy of the simulated electrons, 
        in eV.
        '''
        if 'initial_KE' in kwargs.keys():
            initial_KE = kwargs['initial_KE']
        else:
            initial_KE = 1100
        
        if 'source_diameter' in kwargs.keys():
            source_diameter = kwargs['source_diameter']
        else:
            source_diameter = 300000
            
        if 'source_thickness' in kwargs.keys():
            source_thickness = kwargs['source_thickness']
        else:
            source_thickness = 20
        
        if 'source_angle_dist' in kwargs.keys():
            source_angle_dist = kwargs['source_angle_dist']
        else:
            source_angle_dist = 'Sphere'
                      
        '''Here we define the boundary_shape.'''
        boundary_shape = Disc(sim_radius,2*sample_nozzle_distance)        
        '''Then we instantiate a simulation object and pass it the 
        boundary shape.
        '''
        super().__init__(boundary_shape = boundary_shape, **kwargs)
        
        source_radius = source_diameter / 2  
        source_center = [0,0,- source_thickness / 2]
        
        medium_radius = source_diameter / 2
        medium_thickness = source_diameter / 4
        medium_center = [0,0,- medium_thickness / 2]
        
        ''' This adds a source, i.e. a shape that emits electrons. Here it is a
        Disc. Arguments are disc radius and disc thickness.'''
        self.addSource(source_radius, source_thickness, 
                       center = source_center, 
                       angle_distribution = source_angle_dist) 
        self.initial_KE = initial_KE

        
        '''This defines the scattering medium.'''
        self.addScatteringMedium(Disc(medium_radius, medium_thickness, center = medium_center))
        ''' This sets the parameters of the scatterer.Arguments are:
        denstiy in atoms/nm^3
        inel_factor, and inel_exp which determines the inelastic cross 
        section from a fit to TPP2M and NIST data using 
        sigma = factor / (KE ^ exp)
        el_factor, which determines the elastic cross section
        and atomic number
        One can also pass 'loss_function' as a keyword, where the value
        can either be the filename of a csv file, or a list of lists,
        containing the x-values (energy loss) and y-values (probability)
        of a loss function.
        '''
              
        if 'density' in kwargs.keys():
            density = kwargs['density']
        else:
            density = 10.49
            
        if 'molar_mass' in kwargs.keys():
            molar_mass = kwargs['molar_mass']
        else:
            molar_mass = 47
            
        d = self.convertDensity(density, molar_mass)
        
        
        self.addScatterer(d, 1.551,0.831,0.95,1.12,molar_mass, **kwargs) 
        self.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford', 
                                                        energy = self.initial_KE,
                                                        Z = self.scatterer.Z,
                                                        param = 0.9), 
                               'inelastic':AngleDist(kind = 'Constant')}
        
        '''self.addScatterer(d, 0.004,0,0.00001,0,molar_mass, **kwargs) 
        self.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Constant'), 
                               'inelastic':AngleDist(kind = 'Constant')}'''
        
        
        self.height = 1
        
        self.parameters = {'sim_radius':sim_radius,
                           'sample_nozzle_distance': sample_nozzle_distance,
                           'density':d,
                           'molar_mass':molar_mass,
                           'initial_KE':initial_KE,
                           'source_diameter':source_diameter
                           }
    
    def convertDensity(self, density, mol_mass):
        Av = 6.022E+23
        cm_nm = 1E+21
        d = density / mol_mass * Av / cm_nm
        return d
        

#%% Run the simulation (this could take a while).

if __name__ == '__main__':
    loss_function = 'Ag_loss_fn.csv'
    #loss_function=[[10],[0.5]]
    sim = solidStateSimulation(sample_nozzle_distance = 300000,
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
    
