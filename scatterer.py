# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:52:14 2020

@author: Mark
"""
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
import scipy.stats
from numpy.random import standard_cauchy as cauchy
from shapes import Sphere, Disc
from poisson_distance import Poisson
from generate_angle import AngleDist, Phi
from rotation import length, rotate
from metropolis import Metropolis

#%%

class Scatterer():
    def __init__(self, density, inel_factor, el_factor, el_angle, inel_angle):
        self.inel_factor = inel_factor
        self.el_factor = el_factor
        self.density = density
        # stores the angular probability distributions for the two types of 
        # scattering event
        self.angle_dist = {'elastic':AngleDist(kind = 'Cauchy', width = 1), 
                           'inelastic':AngleDist(kind = 'Constant')}
        #self.angle_dist = {'elastic':AngleDist(kind = 'Constant'), 
        #                   'inelastic':AngleDist(kind = 'Constant')}
        self.avg_loss = 10 # the average amount of kinetic energy lost per
        # inelastic scattering event (in eV)
        self.phi = Phi()
        
    def Scatter(self):
        ''' Determines whether elastic or inelastic scattering occurs.
        Depending on which occurs, it selects the appropriate angular spread
        function. Then it gets a random angle. It returns a list containing
        the kind of scattering event, the distance the electron travelled since
        its last scattering event, and the change in angle upon scattering.
        '''
        d = self.getDistance()
        kind = self.getKind()
        dist = self.angle_dist[kind]
        theta = dist.getAngle()
        phi = self.phi.getAngle()
        return kind, d, theta, phi
        
    def setXSect(self, KE):
        self.inel_xsect = self.inel_factor / KE
        self.el_xsect = self.el_factor / KE
        self.total_xsect = self.inel_xsect + self.el_xsect
        
    def getDistance(self):
        return -1/(self.total_xsect*self.density) * np.log(rand())
    
    def getKind(self):
        if rand() < self.inel_xsect / (self.inel_xsect + self.el_xsect):
            kind = 'inelastic'
        else:
            kind = 'elastic'
        return kind
 