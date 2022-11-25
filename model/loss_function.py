# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:04:05 2020

@author: Mark & Lukas
"""
import matplotlib.pyplot as plt
import numpy as np

try:
    from .probability_distribution import LoadedDistribution, ManualDistribution
    from .samplers import Inversion

except ImportError as e:
    if str(e) == "attempted relative import with no known parent package":
        from probability_distribution import LoadedDistribution, ManualDistribution
        from samplers import Inversion
    else:
        raise

# %%
class LossFunction:
    """The loss function is a probability density function for energy loss."""

    def __init__(self, **kwargs):
        """Construct object."""
        pass

    def getValue(self) -> float:
        """Get a random sample from the distribution."""
        return

class LoadedLossFunction(LossFunction):
    """Loss functions can be loaded from csv. The format col1: x, col2:y."""

    def __init__(self, filename : str):
        """Construct object."""
        self.distribution = LoadedDistribution(filename)
        self.sampler = Inversion(self.distribution)
        self.sampler.current_x = 1
        super().__init__()

    def getValue(self) -> float:
        """Get a random sample from the distribution."""
        return self.sampler.getValue()

class DiscreetLossFunction(LossFunction):
    """The loss function can also be a descreet distribution.

    A list of x values and corresponding y values are manually provided.
    """

    def __init__(self, xy: list):
        """Construct object."""
        self.distribution = ManualDistribution(xy)
        self.sampler = Inversion(self.distribution)
        super().__init__()

    def getValue(self) -> float:
        """Get a random sample from the distribution."""
        return self.sampler.getValue()

#%%
if __name__ == '__main__':
    filename = 'He_loss_fn.csv'
    loss_fn = LoadedLossFunction(filename)
    n = 500000
    values = [loss_fn.getValue() for i in range(n)]
    distribution = np.histogram(values, bins = 150, range=(0,50))
    plt.plot(distribution[1][:-1],distribution[0][:])
