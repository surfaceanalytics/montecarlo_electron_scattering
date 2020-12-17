# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:04:05 2020

@author: Mark
"""

from loaded_distribution import LoadedDistribution, ManualDistribution
from samplers import Metropolis, Inversion
import matplotlib.pyplot as plt
import numpy as np

class LossFunction:
    def __init__(self, **kwargs):
        pass
    
    def getValue(self):
        return          

class LoadedLossFunction(LossFunction):
    """Loss functions can be loaded from csv. The format col1: x, col2:y."""
    
    def __init__(self, filename):
        self.distribution = LoadedDistribution(filename)
        self.sampler = Inversion(self.distribution)
        self.sampler.current_x = 1
        super().__init__()
        
    def getValue(self):
        return self.sampler.getValue()
    
class DiscreetLossFunction(LossFunction):
    def __init__(self, xy):
        self.distribution = ManualDistribution(xy)
        self.sampler = Inversion(self.distribution)
        super().__init__()
        
    def getValue(self):
        return self.sampler.getValue()
    
    
#%%
if __name__ == '__main__':
    filename = 'He_loss_fn_medium.csv'
    loss_fn = LoadedLossFunction(filename)
    n = 500000
    values = [loss_fn.getValue() for i in range(n)]
    distribution = np.histogram(values, bins = 150, range=(0,50))
    plt.plot(distribution[1][:-1],distribution[0][:])
