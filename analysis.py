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
        
        To limit the electrons the the ones that have interesected the
        boundary within a desired acceptance angle, provide the keyword arg
        'accept_angle' = (angle, radius), where angle should be in degrees
        and radius represents the radius of the boundary sphere.
        '''
        max_events = max(self.results[1][:,7])
        if max_events == 0:
            max_events = 1
        hitting_sphere = self.results[1][np.where(self.results[1][:,7]
                                        <max_events)]
        # get a list of the number times an electron was inelastically 
        # scattered for all electrons hitting the sphere
        inel_count = hitting_sphere[:,7]
        # get the pathlengths of all electrons hitting the sphere
        pathlengths = self.results[2][np.where(self.results[1][:,7]
                                        <max_events)]
    
        # This condition restricts the counted electrons to those within a 
        # certain acceptance angle
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
        for i in range(int(max(inel_count))+1):
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
        ''' This plots the total ares under the
        '''
        
        areas = {}
        for k,v in self.profiles.items():
            areas[k] = np.sum(v['counts'])
        self.areas = areas
        if 'show' in args:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
            ax.scatter(list(areas.keys()),list(areas.values()))
            ax.set_xlabel('Number of times inelastically scattered', linespacing=3)
            ax.set_ylabel('Electron count', linespacing=3)
            plt.show()
            
    def plotStartFinish(self, x_lim=800000, n_electrons=500, **kwargs):
        hitting_sphere = self.results[1]
        
        if 'accept_angle' in kwargs.keys():
            acceptance_angle = kwargs['accept_angle'][0] #degrees
            sim_radius = kwargs['accept_angle'][1]
            radius = sim_radius * np.sin(acceptance_angle / 180 * np.pi)
            acceptance = np.array([i for i in hitting_sphere 
                                   if (((i[0]**2 + i[1]**2)
                                   < (radius)**2)
                                    & (i[2] > 0))])
        
        else:
            acceptance = hitting_sphere
        
        if n_electrons > len(acceptance):
            n_electrons = len(acceptance)
        
        fig = plt.figure(figsize=(10,10))
        xlim = x_lim
        ylim = xlim
        zlim = xlim
        ax = fig.gca(projection='3d')
        ax.set_xlabel('x distance', linespacing=3)
        ax.set_ylabel('y distance', linespacing=3)
        ax.set_zlabel('z distance', linespacing=3)
        for l, k in enumerate(acceptance[:n_electrons]):
            V = k
            x, y, z = V[0], V[1], V[2]
    
            ax.scatter3D(x, y, z, c='blue')
        
        for l,k in enumerate(self.results[0][:n_electrons]):
            V = k
            x, y, z = V[0], V[1], V[2]
    
            ax.scatter3D(x, y, z, c='orange')
            
        ax.set_xlim3d(-xlim, xlim)
        ax.set_ylim3d(-ylim,ylim)
        ax.set_zlim3d(-100000,zlim)
        
        ax.set_xticks([-xlim, 0, xlim])
        ax.set_yticks([-ylim, 0, ylim])
        ax.set_zticks([-zlim, 0, zlim])
        
        ax.view_init(10, 90)    
        plt.show()
            


