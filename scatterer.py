# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:52:14 2020

@author: Mark
"""
import numpy as np
from random import random as rand
from random import choice as choice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
from shapes import Sphere, Disc
from generate_angle import AngleDist, Phi

#%%

class Scatterer():
    '''The Scatterer class represents the scattering medium (i.e. some
    atoms, of a particular element, with a specified density).
    The scattering medium can elastically scatter or inelastically scatter
    an electron. Both processes are defined by the elastic and inelastic 
    scattering cross sections.
    The elastic and inelastic scattering processes each have scattering angle 
    distributions. These are probablity distribution functions that generate
    random angles according to the chosen AngleDist class.
    '''
    def __init__(self, density, inel_factor, el_factor, Z):
        self.inel_factor = inel_factor
        self.el_factor = el_factor
        self.density = density
        self.Z = Z
        # stores the angular probability distributions for the two types of 
        # scattering event
        self.angle_dist = {'elastic':AngleDist(kind = 'Cauchy', width = 1), 
                           'inelastic':AngleDist(kind = 'Constant')}
        #self.angle_dist = {'elastic':AngleDist(kind = 'Constant'), 
        #                   'inelastic':AngleDist(kind = 'Constant')}
        self.avg_loss = 10 # the average amount of kinetic energy lost per
        # inelastic scattering event (in eV)
        self.phi = Phi()
        self.loss_fn = [[0.2,10],[0.8,12]] # this represents the loss function
        # it is a list of lists. The elements of the sub-list are [probability
        # energy loss, kinetic energy change in energy loss]
        
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
        ''' Returns a random distance in nm based on an exponential distribution
        from the total scattering cross section (in nm^2)
        '''
        return -1/(self.total_xsect*self.density) * np.log(rand())
    
    def getKind(self):
        if rand() < self.inel_xsect / (self.inel_xsect + self.el_xsect):
            kind = 'inelastic'
        else:
            kind = 'elastic'
        return kind
    
    def getDeltaKE(self):
        c = choice(self.loss_fn)
        r = rand()
        if r < c[0]:
            return c[1]
        else:
            return self.getDeltaKE()
        

