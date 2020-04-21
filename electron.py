# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:47:27 2020

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
from generate_angle import AngleDist
from rotation import length, rotate

#%%

class Electron():
    def __init__(self, source):
        ''' The Electron class represents an electron in the simultation. Its
        main properties are stored inside the attribute called 'vector'.
        Vector is an array to store all the electron's info.
        - the first three elements [0:3] are the x,y,z positions
        - the next three [3:6] are the x, y and z velocities,
        - the next one is time [6]
        - the next one is number of times inelastically scattered [7]
        - the next one is number of time elastically scattered [8]
        The electron needs a Source as an argument, because the electrons need
        to be generated within the volume of the Source.
        '''
        self.vector = np.zeros((9))
        self.source = source
        self.kinetic_energy = 1000 #kinetic energy in eV
        self.setSpeed(self.kinetic_energy) 
        self.theta = AngleDist(kind='Theta')
        self.phi = AngleDist(kind='Phi')
        
    def setSpeed(self, kinetic_energy):
        ''' This function sets the speed of an electron (returned in m/s) using
        kinetic energy as input (units of eV).
        It uses the constants 1.602E-19 Joules per eV and 
        the mass of the electron 9.109E-31 kg.
        '''
        self.initial_speed = np.sqrt(2 * kinetic_energy * 
                                     1.602176E-19 / 9.109383E-31)
        
    def initCoords(self, depth, height):
        ''' This function gets random x,y,z corrdinates from a slice
        of the Source disc. 
        '''
        self.vector[0:3] = self.source.getSlicePosition(depth, height)
        
    def initVelocity(self):
        ''' This  function initializes the electron's velocity. It uses the 
        electron's speed, and generates random polar and azimuthal angles.
        Polar angle is theta, and asimuthal is phi.
        '''
        theta = self.theta.getAngle()
        phi = self.phi.getAngle()
        vx = self.initial_speed * np.sin(theta) * np.cos(phi)
        vy = self.initial_speed * np.sin(theta) * np.sin(phi)
        vz = self.initial_speed * np.cos(theta)
        self.vector[3:6] = np.array([vx,vy,vz])
        
    def changeSpeed(self, delta_KE):
        ''' This fuction is used to change the electron's speed when an
        inelastic scattering event occurs. It takes the argument delta_KE,
        which represents the absolute value of kinetic energy the electron lost
        in the procees (units of eV). It returns a new velocity vector in units
        of m/s.
        '''
        old_v = self.vector[3:6]
        new_KE = self.kinetic_energy - delta_KE
        speed = np.sqrt(2 * new_KE * 
                                     1.602176E-19 / 9.109383E-31)
        new_v = (old_v / np.sqrt(old_v[0]**2 + old_v[1]**2 + old_v[2]**2)
            * speed)
        self.vector[3:6] = new_v
        
#%%

if __name__ == '__main__':
    source = Disc(100,1)
    vals = []
    e = Electron(source)
    for i in range(2000):
        e.initVelocity()
        vals+=[list(e.vector)]
    vals = np.array(vals)
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    xs = vals[:,3]
    ys = vals[:,4]
    zs = vals[:,5]
    ax.scatter(xs, ys, zs, s=20, marker = '.')
  
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.view_init(10, 45)  
    
    plt.show()        
