# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 12:06:39 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
from poisson_distance import Poisson
from probability_distribution import LoadedDistribution

#%%
class Metropolis:
    """A class for using the Metropolis algorithm to return a distribution."""
    
    def __init__(self, distribution):
        self.pdf = distribution
        self.max_x = self.pdf.max_x()
        self.min_x = self.pdf.min_x()
        self.range = self.max_x - self.min_x
        self.n_steps = 10
        self.values = []
        self._initialize()
        
    def _initSteps(self):
        """Determine the step size that should be used."""
        self.step = self.range / self.n_steps
        
    def _initialize(self):
        """Pick an initial x value."""
        self._initSteps()
        self.current_x = rand() * self.range
        self.current_y = self.pdf.y(self.current_x)
        self.values += [self.current_x]
        
    def warmUp(self,n : int):
        """Warm up the Metropolis algorithm."""
        self.run(n)
        self.values = []
    
    def metropolis(self):
        """Pick a value using the Metropolis algorithm."""
        trial_x = self.current_x + (rand() * self.step - self.step / 2)
        if (trial_x < self.min_x) or (trial_x > self.max_x):
            return
        trial_y = self.pdf.y(trial_x)
        R = trial_y / self.current_y
        p = rand()
        if p <= R:
            self.current_y = self.pdf.y(self.current_x)
            self.current_x = trial_x
            self.values += [trial_x]
        else:
            return
            
    def run(self, n: int):
        while len(self.values) < n:
            self.metropolis()
            
    def getValue(self) -> float:
        self.values = []
        self.run(1)
        return self.values[0]
    
#%%
if __name__ == '__main__':
    
    filename = 'He_loss_fn.csv'  
    pdf = LoadedDistribution(filename)
    m = Metropolis(pdf)
    #m = Metropolis(Poisson(0.01,59))
    m.n_steps = 10
    m.current_x = 1
    m.warmUp(2000)

    n_run = 10000
    for i in range(50):
        m.run(n_run)
        v = np.array(m.values)
        dist = 'Energy loss [eV]'
        fig = plt.figure(dpi=400)
        ax = fig.add_subplot(1,1,1,)
        ax.hist(v, bins = 200, color = 'orange')    
    
        plt.xlabel(dist)
        plt.ylabel('Counts')
        plt.show()
        n_run += 10000
        
        dist = 'Energy loss [eV]'
        fig = plt.figure(dpi=400)
        ax = fig.add_subplot(1,1,1,)
        ax.plot(m.pdf.xy[0], m.pdf.xy[1], color = 'blue')    
    
        plt.xlabel(dist)
        plt.ylabel('Relative probability [Arb.U.]')
        plt.show()

        
 
    

