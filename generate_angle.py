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
    ''' a probability dessnity distribution of angles having a Lorentzian
    shape'''
    def __init__(self, **kwargs):
        if 'width' in kwargs.keys():
            self.width = kwargs['width']
        else:
            self.width = 10 / 180 * np.pi
        if 'kind' in kwargs.keys():
            kind = kwargs['kind']
            if kind == 'Cauchy':
                self.distribution = Cauchy(self.width)
            elif kind == 'Gaussian':
                self.distribution = Gaussian(self.width)
            elif kind == 'Spherical':
                self.distribution = Spherical()
            elif kind == 'Phi':
                self.distribution = Phi()
            elif kind == 'Theta':
                self.distribution = Theta()
            elif kind == 'Constant':
                self.distribution = Constant()
            elif kind == 'Rutherford':
                self.distribution = Rutherford()
                    
    def getAngle(self, *args):
        return self.distribution.getAngle(*args)

class Cauchy:
    
    def __init__(self, width):
        self.width = width
        
    def getAngle(self):
        ''' returns a random angle, in degrees, following
        a Lorentzian (Cauchy) distribution'''
        angle = np.abs(cauchy(1)*self.width)
        #if angle > np.pi:
        #    angle = self.getAngle()
        return angle[0]

class Spherical:
    
    def __init__(self):
        pass
    
    def getAngle(self):
        ''' returns a random angle, in degrees, between 0 and 180, following
        a uniform distribution'''
        return rand() * np.pi

class Gaussian:
    
    def __init__(self, width):
        self.width = width
        
    def getAngle(self):
        ''' returns a random angle, in degrees, between 0 and 180, following
        a Gaussian distribution'''
        return np.random.normal(0,self.width)
    
class Phi:
    def __init__(self):
        pass
    
    def getAngle(self):
        return  rand() * 2 * np.pi
    
class Theta:
    def __init__(self):
        pass

    def getAngle(self):
        return np.arccos(1-2*rand())# rand() * np.pi

class Constant:
    def __init__(self):
        pass
    
    def getAngle(self):
        return 0
    
class Rutherford:
    def __init__(self):
        self.alpha = 1000
        pass
    def getAngle(self):
        r = rand()
        angle = np.arccos(1
                        +(2*(1-(1+self.alpha)**r)
                        /self.alpha))
        return angle
        
        
#%%
        
if __name__ == '__main__':
    angle1 = AngleDist(kind = 'Rutherford')
    angle1.distribution.alpha = 3000
    a1 = [angle1.getAngle() for i in range(50000)]
    
    angle2 = AngleDist(kind = 'Rutherford')
    angle2.distribution.alpha = 1000
    a2 = [angle2.getAngle() for i in range(50000)]
    
    angle3 = AngleDist(kind = 'Rutherford')
    angle3.distribution.alpha = 500
    a3 = [angle3.getAngle() for i in range(50000)]
       
    plt.hist(a1, bins = 200, range=(0,5))   
    plt.hist(a2, bins = 200, range=(0,5))   
    plt.hist(a3, bins = 200, range=(0,5))   

