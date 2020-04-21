# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 12:06:39 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
from poisson_distance import Poisson

#%%
class Metropolis:
    def __init__(self, distribution):
        self.pdf = distribution
        self.range = self.pdf.max_x - self.pdf.min_x
        self.n_steps = 15
        self.values = []
        self._initialize()
        
    def _initSteps(self):
        self.step = self.range / self.n_steps
        
    def _initialize(self):
        self.current_x = rand() * self.range
        self.values += [self.current_x]
        
    def warmUp(self,n):
        self._initSteps()
        self.run(n)
        self.values = []
    
    def metropolis(self):
        current_y = self.pdf.y(self.current_x)
        trial_x = self.current_x + (rand() * self.step - self.step / 2)
        trial_y = self.pdf.y(trial_x)
        R = trial_y / current_y
        p = rand()
        if p <= R:
            self.current_x = trial_x
            self.values += [trial_x]
        else:
            pass
            
    def run(self, n):
        self._initSteps()
        while len(self.values) < n:
            self.metropolis()
            
    def getValue(self):
        self.values = []
        self.run(1)
        return self.values[0]
    
#%%
if __name__ == '__main__':
    m = Metropolis(Poisson(0.01,59))
    m.n_steps = 15
    m.warmUp(1000)

    n_run = 50000
    m.run(n_run)
    v = np.array(m.values)
    
    dist = 'Distance'
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,)
    ax.hist(v, bins = 50, color = 'orange')    

    plt.xlabel(dist)
    plt.ylabel('Counts')
    plt.show()