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
    """A class for randomly generating angles with some distribution.
    
    This class contains a probability density function for 
    generating random angles. Several distributions are possible. 
    One choses the distribution using the keyword 'kind'
    Options are: 'Cauchy', 'Gaussian', 'Phi', 'Theta', 'Contant'
    and 'Rutherford'.
    Each distribution has the method 'getAngle()', which returns an angle in
    raidans.
    Some distributions take  attitional arguments. For example, 'Cauchy' and 
    'Gaussian' take 'width' as a keyword argument. 'Rutherford' takes 'energy',
    'Z' and 'param' as keyword arguments.
    """
    
    def __init__(self, **kwargs):
        """Construct the AngleDist object.
        
        Parameters:
        ----------
            **kwargs:
                'width': float
                    An argument defining the angluar spread of the distribution.
                    This argument is only used if the distribution is 'Cauchy' or
                    'Gaussian'.
                'kind': str
                    The kind of angle distribution.
                    One of ['Cauchy', 'Gaussian', 'Phi', 'Theta', 'Constant',
                            'Rutherford', 'Lambert']
        """
        if 'width' in kwargs.keys():
            self.width = kwargs['width']
        else:
            self.width = 10 / 180 * np.pi
        if 'kind' in kwargs.keys():
            kind = kwargs['kind']
            if kind == 'Cauchy':
                self.distribution = Cauchy(**kwargs)
                self.kind = 'Cauchy'
            elif kind == 'Gaussian':
                self.distribution = Gaussian(self.width)
                self.kind = 'Gaussian'
            elif kind == 'Phi':
                self.distribution = Phi()
                self.kind = 'Phi'
            elif kind == 'Theta':
                self.distribution = Theta()
                self.kinf = 'Theta'
            elif kind == 'Constant':
                self.distribution = Constant()
                self.kind = 'Constant'
            elif kind == 'Rutherford':
                self.distribution = Rutherford(**kwargs)
                self.kind = 'Rutherford'
            elif kind == 'Lambert':
                self.distribution = Lambert()
                self.kind = 'Lambert'
                    
    def getAngle(self, *args):
        """Return a random angle, in radians.
                
        Parameters:
        ----------
            *args
        Returns:
        -------
            anlge: float
        """ 
        return self.distribution.getAngle(*args)

class Cauchy:
    """A class fo generating random angles with a Cauchy disribution."""
    
    def __init__(self, **kwargs):
        """Construct a Cauchy object.
        
        Parameters:
        ----------
            **kwargs:
                width: float
                    The width of the distribution.
        """
        if 'width' in kwargs.keys():
            self.width = kwargs['width']
        else:
            self.width = 0.1
        
    def getAngle(self):
        """Return a random angle, in radians.
                
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """ 
        angle = np.abs(cauchy(1)*self.width)
        return angle[0]

class Gaussian:
    """A class fo generating random angles with a Gaussian disribution."""
    
    def __init__(self, width):
        """Construct a Gaussian object.
        
        Parameters:
        ----------
            **kwargs:
                width: float
                    The width of the distribution.
        """
        self.width = width
        
    def getAngle(self):
        """Return a random angle, in radians.
                
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """
        return np.random.normal(0,self.width)
    
class Phi:
    """The Phi class returns an azimuthal angle from a uniform distribution."""
    
    def __init__(self):
        """Construct the object."""
        pass
    
    def getAngle(self):
        """Return a random angle, in radians.
                
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """
        return  rand() * 2 * np.pi
    
class Theta:
    """The Theta class returns a random polar angle.
    
    The distribution generates ploints that are uniform over the surface of a 
    sphere.
    """
    
    def __init__(self):
        """Construct the object."""
        pass

    def getAngle(self):
        """Return a random angle, in radians.
        
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """
        return np.arccos(1-2*rand())

class Constant:
    """The Constant class returns a constant number.
    
    This distributions represent a cases where the angle does not change.
    """
    
    def __init__(self):
        """Construct object."""
        pass
    
    def getAngle(self):
        """Return a constant angle, in radians.
        
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """
        return 0
    
class Rutherford:
    """The Rutherford class returns a random angle from the Rutherford distribution.

    It returns the a random polar angle in radians.
    """
    
    def __init__(self, **kwargs):
        """Construct the object.
        
        Parameters:
        ----------
            **kwargs:
                energy: float
                    The kinetic energy of the scattered electron in eV.
                Z: float
                    The atomic number of the scattering medium
                param: float
                    A parameter that changes the width of the 
                    distribution. Smaller 'param' leads to narrower distributions. 
                    With larger values the distribution converges to the cosin 
                    distribution.
        """
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
            
        self.beta = self._beta(self.energy, self.Z, self.param)

    def getAngle(self):
        """Return a random angle in radians.
        
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """
        r = rand()
        angle = np.arccos((-1-self.beta+r+2*self.beta*r)/(-1-self.beta+r))
        return angle
        
    def _beta(self, energy, Z, param):
        """Generate the beta parameter.
        
        m = 9.10938E-31 # mass of electron in kg
        eV = 1.602176634E-19 # energy of 1 eV in J
        v = np.sqrt(2*eV*energy/m) # velocity of electron in m/s
        a0 = 5.3e-11 # Bohr radius in m
        p = m*v # momentum in kg m/s
        lambda0 = Z**(1/3) / (0.885*a0)
        h = 6.62607E-34 # Planck constant in J s
        beta = 1/4 * (param * h * lambda0 / p)
        
        Parameters:
        ----------
            energy: float
                The kinetic energy of the scattered electron in eV.
            Z: float
                The atomic number of the scatterer
            param: float
                A parameter for determining the width of the distribution
        Returns:
            beta: float
        """
        beta = param * 5.43 * (Z**(2/3))/energy
        return beta
    
class Lambert:
    """The Lambert class returns a random angle from the Lambert distribution.

    It returns the a random polar angle in radians.
    """
    
    def __init__(self, **kwargs):
        """Construct Lambert object."""
        return

    def getAngle(self):
        """Return a random angle in radians.
        
        Parameters:
        ----------
            None
        Returns:
        -------
            anlge: float
        """
        r = rand()
        angle = np.arccos(np.sqrt(1-r))
        return angle
            
#%%
        
if __name__ == '__main__':
    angle1 = AngleDist(kind = 'Rutherford', energy = 1000, Z = 47, param = 10.001)
    a1 = [angle1.getAngle() for i in range(50000)]

    plt.hist(a1, bins = 200, range=(0,5))   

