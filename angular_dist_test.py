# -*- coding: utf-8 -*-
"""
Created on Mon May  4 13:55:05 2020

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
from scatterer import Scatterer
from electron import Electron
import pickle
#%%
def Ruth(Z,E,theta, param):
    beta = param * 5.43 * Z**(2/3) / E
    y = Z**2 / (4 * E**2 * (1-np.cos(theta) + 
                            2 * beta)**2)
    return y
Z = 47
E = 1400
p = 1
ruth = [Ruth(Z,E,i,p) for i in np.arange(0,np.pi,0.01)]
theta = [i for i in np.arange(0,np.pi,0.01)]
a = AngleDist(kind = "Rutherford", energy = 1400, Z=47, param=2)

aa = []

for i in range(200000):
    aa+= [a.getAngle()]


h = np.histogram(aa, bins=100)
r = h[1][:-1]
h1 = [h[0][i]/np.sin(r[i]) for i in range(len(r))]

plt.plot(r[2:],h1[2:]/max(h1[2:]))
plt.plot(theta,ruth/max(ruth))



