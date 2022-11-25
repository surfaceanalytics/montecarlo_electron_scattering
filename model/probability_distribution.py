# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 09:27:58 2020

@author: Lukas
"""
import os

class Distribution():
    """An arbitrary distribution for use in the Metropolis algorithm."""

    def __init__(self):
        """Construct the object."""
        return

    def get_y(self, x: float) -> float:
        """Get the y value given an x value."""
        return

class LoadedDistribution(Distribution):
    """A class that represents a distribution loaded from a file.

    The format of the file should be csv, with the first column as x values
    and the second column as y values.
    """

    def __init__(self, filename: str):
        """Construct object."""
        loader = CSVLoader()
        self.xy = loader.Load(filename)
        super().__init__()
        self._construct_CDF()

    def get_y(self,x:float) -> float:
        """Return a sample from the probability distribution."""
        _x = self.xy[0]
        '''Get the index for the value closest to x.'''
        idx = min(range(len(_x)), key=lambda i: abs(_x[i]-x))
        y = self.xy[1][idx]
        return y

    def max_x(self):
        """Return the maximum x-value in the distribution."""
        return max(self.xy[0])

    def min_x(self):
        """Return the minimum x-value in the distribution."""
        return min(self.xy[0])

    def _construct_CDF(self):
        y_sum = 0
        for y in self.xy[1]:
            y_sum += y
        self.xy[1] = [y/y_sum for y in self.xy[1]]

class ManualDistribution(Distribution):
    """A class that represents a distribution loaded from a file.

    The format of the file should be csv, with the first column as x values
    and the second column as y values.
    """

    def __init__(self, xy: list):
        """Construct object.

        The argument xy should be a list of lists, with the first item being
        a list of x values and the second item being a list of y values.
        """
        self.xy = xy
        super().__init__()
        self._construct_CDF()

    def get_y(self, x:float) -> float:
        """Return a sample from the probability distribution."""
        _x = self.xy[0]
        '''Get the index for the value closest to x.'''
        idx = min(range(len(_x)), key=lambda i: abs(_x[i]-x))
        y = self.xy[1][idx]
        return y

    def max_x(self):
        """Return the maximum x-value in the distribution."""
        return max(self.xy[0])

    def min_x(self):
        """Return the minimum x-value in the distribution."""
        return min(self.xy[0])

    def _construct_CDF(self):
        """Construct the cumulative distribution function."""
        y_sum = 0
        for y in self.xy[1]:
            y_sum += y
        self.xy[1] = [y/y_sum for y in self.xy[1]]


class CSVLoader():
    """A loader class for csv files of x and y values.

    Format should be first column: x-values, second column: y-value.
    """

    def __init__(self):
        """Construct object."""
        pass

    def Load(self, filename: str) -> list:
        """Load the csv file into a list of lists."""
        x=[]
        y=[]
        parentpath = os.path.join(
            os.getcwd().split("montecarlo")[0],
            "montecarlo"
            )
        filepath = os.path.join(*[parentpath,"data",filename])
        with open(filepath, 'r') as f:
            for line in f.readlines():
                try:
                    x += [float(line.split(',')[0])]
                    y += [float(line.split(',')[1])]
                except:
                    continue
        xy = [x,y]
        return xy

if __name__=='__main__':
    filename = 'He_loss_fn.csv'

    dist = LoadedDistribution(filename)

    print(dist.get_y(50))


