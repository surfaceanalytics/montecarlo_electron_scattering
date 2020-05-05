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
#%%

class Analysis:
    def __init__(self, results):
        ''' The analysis class takes a complete results set in as 
        '''
        self.results = results
        
    def sphere_intersections(self):
        # first pick only the electrons that hit the upper half of the sphere
        max_events = max(self.results[1][:,7])
        hitting_sphere = self.results[np.where(self.results[1][:,7]<max_events)]
        # then get the final position and velocity of these electrons
        return hitting_sphere
    
    def pathHistogram(self,*args, bins=100, **kwargs):
        '''Gets the max number fo times inelastically scattered
        this should represent the cutoff for number of times scattering 
        in the simublation before the path calculation was terminated
        
        Use the argument 'show' to show the plot after, and the keyword
        argument 'x_limits' with a tuple (x_min, x_max) to adjust the axis
        limits
        
        To limite the electrons the the ones that have interesected the
        boundary within a desired acceptance angle, provide the keayword arg
        'accept_angle' = (angle, radius), where angle should be in degrees
        and radius represents the radius of the boundary sphere.
        '''
        max_events = max(self.results[1][:,7])
        hitting_sphere = self.results[1][np.where(self.results[1][:,7]
                                        <max_events)]
        # get a list of the number times an electron was inelastically 
        # scattered for all electrons hitting the sphere
        inel_count = hitting_sphere[:,7]
        # get the pathlengths of all electrons hitting the sphere
        pathlengths = self.results[2][np.where(self.results[1][:,7]
                                        <max_events)]
    
        if 'accept_angle' in kwargs.keys():
            acceptance_angle = kwargs['accept_angle'][0] #degrees
            radius = kwargs['accept_angle'][1] * np.sin(acceptance_angle / 
                           180 * np.pi)
            hitting_sphere = np.array([i for i in hitting_sphere if (((i[0]**2
                                       +i[1]**2) < (radius)**2) & (i[2] > 0))])
            inel_count = hitting_sphere[:,7]
    
            pathlengths = np.array([self.results[2][i] 
                    for i in range(len(self.results[2]))                        
                    if ((self.results[1][i][0]**2
                            + self.results[1][i][1]**2
                            < (radius)**2)
                            & (self.results[1][i][2] > 0))])

        profiles = {}
        for i in range(int(max(inel_count))):
            hist = np.histogram(pathlengths[np.where(inel_count == i)], 
                                            bins=100)
            profiles[i] = {'counts':hist[0], 'pathlength':hist[1][:-1]}
        
        if 'show' in args:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
            for k,v in profiles.items():
                ax.plot(v['pathlength'], v['counts'])
            if 'x_limits' in kwargs.keys():
                limits = kwargs['x_limits']
                ax.set_xlim(limits)
            ax.set_xlabel('Pathlength [nm]', linespacing=3)
            ax.set_ylabel('Electron count', linespacing=3)
            plt.xticks(rotation=90)
            plt.show()
        self.profiles = profiles
    
    def areasUnderProfiles(self, *args):
        areas = {}
        for k,v in self.profiles.items():
            areas[k] = np.sum(v['counts'])
        self.areas = areas
        if 'show' in args:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
            ax.plot(list(areas.keys()),list(areas.values()))
            ax.set_xlabel('Number of times inelastically scattered', linespacing=3)
            ax.set_ylabel('Electron count', linespacing=3)
            plt.show()
            
        
            
#%%
if __name__ == '__main__':
    sim = Simulation(800000) 
    # argument is the radius of the spherical simulation evnrionment
    sim.addSource(100000,5000) # source is a shape that emits electrons
    # arguments are radius and thickness (in nm) of a disc
    sim.addScatteringMedium(Sphere(800000))
    sim.addScatterer(9.66E-5,0.527,0.85,0.5,1.01,4) 
    # arguments are denstiy in atoms/nm^3
    # inel_factor, and inel_exp which determines the inelastic cross section
    # from a fit to TPP2M and NIST data using sigma = facter / (KE * exp)
    # el_factor, which determines the elastic cross section
    # and atomic number
    sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford', 
                                                    energy = sim.initial_KE,
                                                    Z = sim.scatterer.Z,
                                                    param = 0.06), 
                           'inelastic':AngleDist(kind = 'Constant')}
    #set the high of the source slice
    sim.height = 1

#%% 

# potentially very long calculation
if __name__ == '__main__':

    sim.simulateMany(100000, 'start finish')
    results = sim.start_finish

#%%
if __name__ == '__main__':
    A = Analysis(results)
    A.pathHistogram('show',bins=25)        
    A.areasUnderProfiles('show')

    
        
    #%%
    # This plots the average number of times being inelastically scattered versus
    # depth
    # first pick only the electrons that hit the upper half of the sphere
    hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]
    # then get the number of times these electrons were inelastically scattered
    scatter_count = hitting_sphere[:,-1,7]
    # then get the initial position of these electrons
    first_position = hitting_sphere[:,0,:] 
    
    # change scatter count column 
    first_position[:,7] = scatter_count
    
    # find min and max z positions
    min_z = np.min(first_position[:,2])
    max_z = np.max(first_position[:,2])
    step_size = 0.5
    steps = np.arange(-50,0,step_size)
    avg_depth = []
    cts_depth = []
    
    for s in steps:
        avg = np.mean(first_position[np.where((first_position[:,2] < s + step_size)
        &(first_position[:,2] >= s)),7])
        cts = len((first_position[np.where((first_position[:,2] < s + step_size)
        &(first_position[:,2] >= s)),7][0,:]))
        avg_depth += [[s,avg]]
        cts_depth += [[s,cts]]
    avg_depth = np.array(avg_depth)
    cts_depth = np.array(cts_depth)
    
    # this figure shows the number of times inelastically scattered versus the
    # average depth below the surface for electrons being scattered that many times
    fig = plt.figure(figsize=(5,5))
    ax = fig.gca()
    ax.plot(avg_depth[:,0], avg_depth[:,1])
    ax.set_xlabel('Depth [nm]', linespacing=3)
    ax.set_ylabel('Number of times inelastically scattered', linespacing=3)
    plt.show()
    
    # This plot shows the number of electrons reaching the boundary versus depth
    fig = plt.figure(figsize=(5,5))
    ax = fig.gca()
    ax.plot(cts_depth[:,0], cts_depth[:,1])
    ax.set_xlabel('Depth [nm]', linespacing=3)
    ax.set_ylabel('Number of electrons reaching boundary', linespacing=3)
    plt.show()
    
    
    
    #%%
    # this shows the number of electrons reachng the boundary versus depth for 
    # a given number of scattering events
    # the plots should have a Poisson distribution
    
    # first pick only the electrons that hit the upper half of the sphere
    hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]
    # then get the number of times these electrons were inelastically scattered
    scatter_count = hitting_sphere[:,-1,7]
    # then get the initial position of these electrons
    first_position = hitting_sphere[:,0,:] 
    
    # change scatter count 
    first_position[:,7] = scatter_count
    
    # This plots shows the number of electrons reaching the boundary for a given 
    # n. The plots should follow the Poisson distribution
    fig = plt.figure(figsize=(5,5))
    ax = fig.gca()
    
    all_distributions = []
    for n in range(12):
        single_n = first_position[np.where(first_position[:,7] == n)]
        
        # find min and max z positions
        min_z = np.min(first_position[:,2])
        max_z = np.max(first_position[:,2])
        step_size = 0.5
        steps = np.arange(-40,0,step_size)
        cts_depth = []
        
        for s in steps:
            cts = len((single_n[np.where((single_n[:,2] < s + step_size)
            &(single_n[:,2] >= s)),7][0,:]))
            cts_depth += [[s,cts]]
        cts_depth = np.array(cts_depth)
        all_distributions +=[cts_depth[:,1]]
    
        
        ax.plot(cts_depth[:,0], cts_depth[:,1])
    
    ax.set_xlabel('Depth [nm]', linespacing=3)
    ax.set_ylabel('Number of electrons reaching boundary', linespacing=3)
    plt.show()
    
    all_distributions = np.array(all_distributions).T
    
    sums = np.sum(all_distributions, axis=0)