# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 13:57:10 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import standard_cauchy as cauchy
from numpy.random import random as rand

#%%
class AngleDist():
    ''' This class represents a probability density function for choices of 
    random angles. Several distributions are possible. One choses the 
    distribution using the keyword 'kind'
    Options are: 'Cauchy', 'Gaussian', 'Phi', 'Theta', 'Contant'
    and 'Rutherford'
    Each distribution has the method 'getAngle()', which returns an angle in
    raidans.
    Some distributions take  attitional arguments. For example, 'Cauchy' and 
    'Gaussian' take 'width' as a keyword argument. 'Rutherford' takes 'energy',
    'Z' and 'param' as keyword arguments.
    '''
    def __init__(self, **kwargs):
        if 'width' in kwargs.keys():
            self.width = kwargs['width']
        else:
            self.width = 10 / 180 * np.pi
        if 'kind' in kwargs.keys():
            kind = kwargs['kind']
            if kind == 'Cauchy':
                self.distribution = Cauchy(**kwargs)
            elif kind == 'Gaussian':
                self.distribution = Gaussian(self.width)
            elif kind == 'Phi':
                self.distribution = Phi()
            elif kind == 'Theta':
                self.distribution = Theta()
            elif kind == 'Constant':
                self.distribution = Constant()
            elif kind == 'Rutherford':
                self.distribution = Rutherford(**kwargs)
                    
    def getAngle(self, *args):
        return self.distribution.getAngle(*args)

class Cauchy:
    
    def __init__(self, **kwargs):
        if 'width' in kwargs.keys():
            self.width = kwargs['width']
        else:
            self.width = 0.1
        
    def getAngle(self):
        ''' returns a random angle, in degrees, following
        a Lorentzian (Cauchy) distribution'''
        angle = np.abs(cauchy(1)*self.width)
        #if angle > np.pi:
        #    angle = self.getAngle()
        return angle[0]

class Gaussian:
    def __init__(self, width):
        self.width = width
        
    def getAngle(self):
        return np.random.normal(0,self.width)
    
class Phi:
    ''' This class returns an azimuthal angle from a uniform distribution.
    '''
    def __init__(self):
        pass
    
    def getAngle(self):
        return  rand() * 2 * np.pi
    
class Theta:
    ''' This class represents a polar angle using a distribution function
    that is uniform over the surface of a sphere.
    '''
    def __init__(self):
        pass

    def getAngle(self):
        return np.arccos(1-2*rand())# rand() * np.pi

class Constant:
    ''' This clas srepresents an angular distribution, where the angle does 
    not change
    '''
    def __init__(self):
        pass
    
    def getAngle(self):
        return 0
    
class Rutherford:
    ''' This class represents a probability density function for the Rutherford
    scattering profile.
    It returns the a random polar angle in radians.
    It takes three keyword arguments: 
        1) 'energy' represents the electron's kinetic energyy in eV
        2) 'Z', represents the atomic number of the scattering medium
        3) 'param' represents a parameter that changes the width of the 
        distribution. Smaller 'param' leads to narrower distributions. With 
        larger 'param' the distribution converges to the cosin distribution.
    '''
    def __init__(self, **kwargs):
        if 'energy' in kwargs.keys():
            self.energy = kwargs['energy']
        else: 
            self.energy = 1000
        if 'Z' in kwargs.keys():
            self.Z = kwargs['Z']
        else:
            self.Z = 1
        if 'param' in kwargs.keys():
            self.param = kwargs['param']
        else:
            self.param = 0.1
            
        self.beta = self.beta(self.energy, self.Z, self.param)

    def getAngle(self):
        r = rand()
        angle = np.arccos(1-(2*self.beta*r/(1+self.beta - r)))
        return angle
        
    def beta(self, energy, Z, param):
        m = 9.10938E-31 # mass of electron in kg
        eV = 1.602176634E-19 # energy of 1 eV in J
        v = np.sqrt(2*eV*energy/m) # velocity of electron in m/s
        a0 = 5.3e-11 # Bohr radius in m
        p = m*v # momentum in kg m/s
        lambda0 = Z**(1/3) / (0.885*a0)
        h = 6.62607E-34 # Planck constant in J s
        beta = 1/4 * (param * h * lambda0 / p)
        return beta
            
#%%
        
if __name__ == '__main__':
    angle1 = AngleDist(kind = 'Rutherford', energy = 1000, Z = 47, param = 10.001)
    a1 = [angle1.getAngle() for i in range(50000)]

    plt.hist(a1, bins = 200, range=(0,5))   

