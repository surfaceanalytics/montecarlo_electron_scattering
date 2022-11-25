# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:55:25 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
import time

try:
    from .probability_distribution import LoadedDistribution, ManualDistribution

except ImportError as e:
    if str(e) == "attempted relative import with no known parent package":
        from probability_distribution import LoadedDistribution, ManualDistribution
    else:
        raise

#%%

class Sampler:
    """The Sampler is a class for different ways of sampling distributions."""

    def __init__(self, distribution):
        """Construct object.

        Parmeters
        ---------
            distribution: list
                A nested list of lists. The upper list contains two items.
                These are a list of x values and a list of y values.
        """
        self.distribution = distribution

class Inversion(Sampler):
    """The is an inversion sampler for an arbitrary function of xy pairs.

    the X-values represent some observable, and the y values represent the
    probability of measuring that observable.
    This class first normalizes the y-values, so that the total of y-values
    sums to 1. Then it builds a cumulative distribution function.
    Then it uses a uniform random number between 0-1. Gives that as the
    cumulative probability value, gets the index of the value in return,
    and uses the index to retrieve the corresponding observable value.
    """

    def __init__(self, distribution):
        """Construct the object."""
        super().__init__(distribution)
        self._buildCDF()

    def _buildCDF(self):
        self.cumulative_dist = []
        for y in self.distribution.xy[1]:
            if len(self.cumulative_dist) != 0:
                self.cumulative_dist += [self.cumulative_dist[-1] + y]
            else:
                self.cumulative_dist = [y]

    def getValue(self) -> float:
        """Randomly select observable value."""
        r = rand()
        _x = self.cumulative_dist
        '''Get the index for the value closest to x.'''
        idx = min(range(len(_x)), key=lambda i: abs(_x[i]-r))
        return self.distribution.xy[0][idx]

class Metropolis(Sampler):
    """A class for using the Metropolis algorithm to return a distribution."""

    def __init__(self, distribution):
        """Construct object."""
        super().__init__(distribution)
        self.max_x = self.distribution.max_x()
        self.min_x = self.distribution.min_x()
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
        self.current_y = self.distribution.get_y(self.current_x)
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
        trial_y = self.distribution.get_y(trial_x)
        R = trial_y / self.current_y
        p = rand()
        if p <= R:
            self.current_y = self.distribution.get_y(self.current_x)
            self.current_x = trial_x
            self.values += [trial_x]
        else:
            return

    def run(self, n: int):
        """Run the algorithm n times."""
        while len(self.values) < n:
            self.metropolis()

    def getValue(self) -> float:
        """Get one value."""
        self.values = []
        self.run(1)
        return self.values[0]

#%% Usage test
if __name__ == '__main__':

    # test the Metropolis sampler
    filename = 'He_loss_fn.csv'
    pdf = LoadedDistribution(filename)
    m = Metropolis(pdf)
    #m = Metropolis(Poisson(0.01,59))
    m.n_steps = 10
    m.current_x = 1
    m.warmUp(2000)

    start_time = time.time()
    m.run(10000)
    stop_time = time.time()
    print(stop_time - start_time)

    v = np.array(m.values)
    dist = 'Energy loss [eV]'
    fig = plt.figure(dpi=400)
    ax = fig.add_subplot(1,1,1,)
    ax.hist(v, bins = 200, color = 'orange')

    plt.xlabel(dist)
    plt.ylabel('Counts')
    plt.show()

    # Test the inversion sampler
    filename = 'He_loss_fn.csv'
    pdf = LoadedDistribution(filename)
    inv = Inversion(pdf)

    start_time = time.time()
    values = []
    for i in range(10000):
        values += [inv.getValue()]
    stop_time = time.time()
    print(stop_time - start_time)

    dist = 'Energy loss [eV]'
    fig = plt.figure(dpi=400)
    ax = fig.add_subplot(1,1,1,)
    ax.hist(values, bins = 200, color = 'orange')

    plt.xlabel(dist)
    plt.ylabel('Counts')
    plt.show()