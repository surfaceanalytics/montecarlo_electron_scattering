# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 10:13:35 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt

#%%
class Poisson():
    def __init__(self,sigma,rho):
        self.n = 1
        self.sigma = sigma
        self.rho = rho
        self.cut_off = 1E-9
        self._max_y()
        self._set_nmfp()
        self.min_x = 0
        self.max_x = self.n_mfp * self.mfp
        self.fractionReject()

    def y(self,x):
        ''' returns the probability desity function of a particle being 
        scattered 1 time over a distance of x'''
        y = ((1/np.math.factorial(self.n)) 
            * ((x * self.rho * self.sigma)**self.n) 
            * np.exp(-1 * x * self.rho * self.sigma))
        return y
    
    def min_x(self):
        return self.min_x
    
    def max_x(self):
        return self.max_x
    
    def max_y(self):
        self._max_y()
        return self.max_y
    
    def cdf(self, x):
        ''' cumulative distributino function of the Poisson distribution'''
        a = self.rho * self.sigma
        cdf = 1/a - (a*x+1)*np.exp(-1*a*x)/a 
        return cdf

    def _max_y(self):
        ''' finds the maximum value of y (probility)
        '''
        self.max_y = self.y(1/(self.sigma*self.rho))

    def _set_nmfp(self):
        ''' determines the number of mfp's to use in order to reach the desired
        cutoff. The cutoff is a probability, such that distance that have a 
        lower probability than our cutoff are not sampled.
        '''
        prob = 1
        mfp = 1/(self.sigma * self.rho)
        i = 1
        while prob > self.cut_off:
            x = i * mfp
            prob = self.y(x)
            n_mfp = i
            i += 1
        self.n_mfp = n_mfp
        self.mfp = mfp
        
    def fractionReject(self):
        ''' This determines the fraction of attempts of seleceting from the 
        probability distribution that weill be successful
        '''
        self._max_y()
        self._set_nmfp()
        rect_area = self.n_mfp * self.mfp* self.max_y    
        poiss_area = self.cdf(self.n_mfp*self.mfp)    
        self.ratio = poiss_area / rect_area
        
    def _rand(self):
        ''' this picks a random x and y value, within the range of (0, max(x))
        and (0, max(y))
        '''
        rand_x = rand() * (self.n_mfp * self.mfp)
        rand_y = rand() * self.max_y
        return rand_x, rand_y
    
    def _compare(self):
        ''' this compares the randomly selected x and y values with the Poisson
        distribution. If the value of the Poisson distribution a the random
        value of x is less than the value of y, then the x value is successful.
        '''
        x, y = self._rand()
        if self.y(x) >= y:
            return True, x
        else:
            return False, None

    def getDistance(self):
        ''' This takes the successful value of x and returns it along with the
        number of attempts needed to get it.
        '''
        accept = False
        attempts = 0
        while not accept:
            attempts +=1
            accept, x = self._compare()
        return x, attempts

#%%
if __name__ == '__main__':        
    p = Poisson(0.01,59)
    
    trials = []
    for i in range(500000):
        x,a = p.getDistance()
        trials += [[x,a]]
        
    avg_attempts = np.average(np.array(trials)[:,1])
    
    distances = np.array(trials)[:,0]
    '''
    p = Poisson(0.022,59)
    
    trials = []
    for i in range(100000):
        x,a = p.getDistance()
        trials += [[x,a]]
        
    avg_attempts = np.average(np.array(trials)[:,1])
    
    distances1 = np.array(trials)[:,0]'''
    
    '''
    p = Poisson(0.032,59)
    
    trials = {'el':[],'inel':[]}
    for i in range(100000):
        x,a = p.getDistance()
        r = rand()
        if r < (0.022 / (0.01+0.022)):
            trials['el'] += [[x,a]]
        else:
            trials['inel'] += [[x,a]]

    el_distances = (np.array(trials['el'])[:,1])
    inel_distances = (np.array(trials['inel'])[:,1])
    '''
    
    
    dist = 'Distance'
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,)
    ax.hist(distances, bins = 50, color = 'blue')


    plt.xlabel(dist)
    plt.ylabel('Counts')
    plt.show()

