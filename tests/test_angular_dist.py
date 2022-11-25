# -*- coding: utf-8 -*-
"""
Created on Mon May  4 13:55:05 2020

@author: Mark
"""
import os
import pathlib
os.chdir(pathlib.Path(__file__).parent.parent.parent)

import numpy as np
import matplotlib.pyplot as plt

from montecarlo.model.generate_angle import AngleDist
#%%
def Ruth(Z,E,theta, param):
    beta = param * 5.43 * Z**(2/3) / E
    y = Z**2 / (4 * E**2 * (1-np.cos(theta) +
                            2 * beta)**2)
    return y
Z = 47
E = 400
p = 0.9
ruth = [Ruth(Z,E,i,p) for i in np.arange(0,np.pi,0.01)]
theta = [i for i in np.arange(0,np.pi,0.01)]
a = AngleDist(kind = "Rutherford", energy = 400, Z=47, param=0.9)

aa = []

for i in range(10000000):
    aa+= [a.getAngle()]


h = np.histogram(aa, bins=100)
r = h[1][:-1]
h1 = [h[0][i]/np.sin(r[i]) for i in range(len(r))]

plt.plot(r[1:],h1[1:]/max(h1[1:]))
plt.plot(theta,ruth/max(ruth))



