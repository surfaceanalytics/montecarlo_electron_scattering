# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:23:29 2020

@author: Mark
"""
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
import scipy.stats
from numpy.random import standard_cauchy as cauchy
#%%

class Sphere():
    def __init__(self,r):
        self.r = r
        
    def xy(self, z):
        xyz = []
        if z**2 <= self.r**2:
            for theta in np.arange(-np.pi,np.pi,np.pi*2/16):
                R = np.sqrt((self.r)**2 - z**2)
                x = R * np.cos(theta)
                y = R * np.sin(theta)
                xyz += [[x,y,z]]
        else:
            print('z not within sphere')
        return xyz
    
    def inside(self, x,y,z):
        inside = True
        if np.sqrt(x**2 + y**2 + z**2) > self.r:
            inside = False
        else:
            inside = True
        return inside
    
    def getIntersection(self,vector):
        x = vector[0]
        y = vector[1]
        z = vector[2]
        x_norm = x / (np.sqrt(x**2 + y**2 + z**2))
        y_norm = y / (np.sqrt(x**2 + y**2 + z**2))
        z_norm = z / (np.sqrt(x**2 + y**2 + z**2))
        x_inter = x_norm * self.r
        y_inter = y_norm * self.r
        z_inter = z_norm * self.r
        intersection = np.array([x_inter,y_inter,z_inter])
        return intersection
            
class Disc():
    def __init__(self, r, h):
        self.r = r
        self.h = h
        
    def inside(self, position):
        x = position[0]
        y = position[1]
        z = position[2]
        inside = True
        if ((z < -self.h)
        | (z > 0)
        | (np.sqrt(x**2  + y**2) > self.r**2)):
            inside = False
        else:
            inside = True
        return inside
    
    def getRandPosition(self):
        theta = rand()*2*np.pi - np.pi
        r = rand() * self.r
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = rand()*self.h - self.h
        return np.array([x,y,z])
    
    def getSlicePosition(self, depth, height):
        ''' This function generates a random position iside of a 'slice'
        The slice represents a part of the Source disc. It has the same radius
        as the Source, but it has its own height and depth below the top 
        surface of the disc.
        '''
        theta = rand()*2*np.pi - np.pi
        r = rand() * self.r
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = rand()*(-height) - depth
        return np.array([x,y,z])

#%%
        
if __name__ == '__main__':
    d = Disc(5000,100)  
    d_xyz = []
    for i in range(10000):
        d_xyz += [d.getSlicePosition(0.0001,0.01)]
    d_xyz = np.array(d_xyz) 
    d_xyz1 = []
    for i in range(10000):
        d_xyz1 += [d.getSlicePosition(99.9,0.01)]
    d_xyz1 = np.array(d_xyz1)  

    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111, projection='3d')
    xs = d_xyz[:,0]
    ys = d_xyz[:,1]
    zs = d_xyz[:,2]
    ax.scatter(xs, ys, zs)
    
    xs = d_xyz1[:,0]
    ys = d_xyz1[:,1]
    zs = d_xyz1[:,2]
    ax.scatter(xs, ys, zs)
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
     
    s = Sphere(5000)
    s_xyz = []
    for z in np.arange(-s.r,s.r,s.r/4):
        s_xyz += s.xy(z)
        
    s_xyz = np.array(s_xyz)    
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    xs = s_xyz[:,0]
    ys = s_xyz[:,1]
    zs = s_xyz[:,2]
    ax.scatter(xs, ys, zs, marker='.')
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show()
    
    
    #%%
    xyz_int = []
    for i in range(2000):
        xyz_int += [s.getIntersection([rand()-0.5,rand()-0.5,rand()-0.5])]
    xyz_int = np.array(xyz_int)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs = xyz_int[:,0]
    ys = xyz_int[:,1]
    zs = xyz_int[:,2]
    
    ax.scatter(xs,ys,zs, color='orange', marker='.', s=20)
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show()
    
    #%%
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs = d_xyz[:,0]
    ys = d_xyz[:,1]
    zs = d_xyz[:,2]
    ax.scatter(xs, ys, zs, s=20, marker = '.')
    xs1 = s_xyz[:,0]
    ys1 = s_xyz[:,1]
    zs1 = s_xyz[:,2]
    ax.scatter(xs1,ys1,zs1, color='orange', marker='.', s=20)
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show()
    