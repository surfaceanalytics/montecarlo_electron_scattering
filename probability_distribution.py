# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 09:27:58 2020

@author: Mark
"""

class Distribution():
    """An arbitrary distribution for use in the Metropolis algorithm."""
    
    def __init__(self):
        """Construct the object."""
        return
    
    def get_y(self, x: float) -> float:
        """Get the y value given an x value."""
        return
    
class LoadedDistribution(Distribution):
    """ A class that represents a distribution loaded from a file.
    
    The format of the file should be csv, with the first column as x values
    and the second column as y values.
    """
    
    def __init__(self, filename: str):
        """Construct object."""
        loader = CSVLoader()
        self.xy = loader.Load(filename)  
        super().__init__()
        self.normalize_y()
        
    def get_y(self,x:float) -> float:
        _x = self.xy[0]
        '''Get the index for the value closest to x.'''
        idx = min(range(len(_x)), key=lambda i: abs(_x[i]-x))
        y = self.xy[1][idx]
        return y
        
    def max_x(self):
        return max(self.xy[0])
    
    def min_x(self):
        return min(self.xy[0])
    
    def normalize_y(self):
        y_sum = 0
        for y in self.xy[1]:
            y_sum += y
        self.xy[1] = [y/y_sum for y in self.xy[1]]

class ManualDistribution(Distribution):
    """ A class that represents a distribution loaded from a file.
    
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
        self.normalize_y()
        
    def get_y(self, x:float) -> float:
        _x = self.xy[0]
        '''Get the index for the value closest to x.'''
        idx = min(range(len(_x)), key=lambda i: abs(_x[i]-x))
        y = self.xy[1][idx]
        return y
               
    def max_x(self):
        return max(self.xy[0])
    
    def min_x(self):
        return min(self.xy[0])
    
    def normalize_y(self):
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
        with open(filename, 'r') as f:
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
    
    print(dist.y(50))

    