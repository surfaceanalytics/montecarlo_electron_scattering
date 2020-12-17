# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 16:12:25 2020

@author: Mark
"""
from shapes import Sphere, Disc
from generate_angle import AngleDist

from numpy.random import rand
import numpy as np

class Source(Disc):
    def __init__(self, r, h, **kwargs):
        super().__init__(r,h,**kwargs)
        
        if 'angle_distribution' in kwargs.keys():
            if kwargs['angle_distribution'] == 'Sphere':
                self.theta = AngleDist(kind='Theta')
                self.phi = AngleDist(kind='Phi')
            elif kwargs['angle_distribution'] == 'Lambert':
                self.theta = AngleDist(kind='Lambert')
                self.phi = AngleDist(kind='Phi')
            elif kwargs['angle_distribution'] == 'Constant':
                self.theta = AngleDist(kind='Constant')
                self.phi = AngleDist(kind='Phi')
        else:
            self.theta = AngleDist(kind='Lambert')
            self.phi = AngleDist(kind='Phi')
            
        if 'initial_KE' in kwargs.keys():
            self.initial_KE = kwargs['initial_KE']
        else:
            self.initial_KE = 1100
            
        if 'emission_line' in kwargs.keys():
            self.emission_line =kwargs['emission_line']
        else:
            self.emission_line = ('d', 5.98)
            
        self.E_distribution = EnergyDistribution(self.initial_KE, 
                                                 self.emission_line[0],
                                                 split = self.emission_line[1])
        
    def getKE(self):
        return self.E_distribution.getKE()
        
class EnergyDistribution:
    
    def __init__(self, initial, line_type, **kwargs):
        self.width = 0.3
        if 'split' in kwargs.keys():
            split = kwargs['split']
        else:
            split = 5.98
            
        if line_type == 's':    
            self.E = [[1],[initial]]
        elif line_type == 'p':
            self.E = [[2/3,1/3],[initial,initial-split]]
        elif line_type == 'd':
            self.E = [[3/5,2/5],[initial, initial-split]]
        elif line_type == 'f':
            self.E = [[4/7, 3/7],[initial, initial-split]]
            
    def getKE(self):
        r = rand()
        if r <= self.E[0][0]:
            KE = self.E[1][0]
        else:
            KE = self.E[1][1]
        return np.random.normal(KE, self.width)