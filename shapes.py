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
#%%

class Sphere():
    """ Sphere is a class that represents the simulation environment. Electrons
    are generated form a 'source' shape inside of the sphere, and eventually 
    intersect the sphere during the simulation.
    The Sphere class has a method to generate x and y positions, given a z 
    It has a method that checks if a point[x,y,z] is inside the sphere. This 
    method returns a boolean.
    It also has a method to return the intersection of a vector with the sphere.
    The method takes a vector, that represents direction, as input, and returns
    x,y,z as outputs, representing the intersection with the sphere.    
    """
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
    
    def inside(self, position):
        ''' This checks if a given point is inside of the shape
        '''
        x = position[0]
        y = position[1]
        z = position[2]
        inside = True
        if np.sqrt(x**2 + y**2 + z**2) > self.r:
            inside = False
        else:
            inside = True
        return inside
    
    def getIntersection(self,vector):
        ''' This finds the intersection point on the shape's surface for a 
        vectors whose start point is at the center of the shape.
        '''
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
    
    def findIntersect(self, position, direction):
        ''' This function finds the intersection between a line and a shape's
        surface. The line has a psition and a direction.The function returns
        the coordinates of the intersection.
        '''
        v = direction / np.linalg.norm(direction)
        p = position
        dot = np.dot(v,p)
        d1 = -dot + np.sqrt(dot**2-(np.linalg.norm(p)**2-self.r**2))
        d2 = -dot - np.sqrt(dot**2-(np.linalg.norm(p)**2-self.r**2))
        d = np.max([d1,d2]) 
        intersect = p + (v*d)
        return intersect
            
class Disc():
    ''' The Disc class is used for the shape where electrons are generated and
    where they are scattered. It is defined by the attributes r (radius in nm)
    and h (height in nm).
    It has a method to return a boolean, determining if a point is inside the
    volums of the disc or not.
    It has a method to return a random position from inside the volume of the 
    disc.
    It has a method that generates a random position from a slice inside of the
    disc.
    '''
    def __init__(self, r, h, **kwargs):
        self.r = r
        self.h = h
        if 'center' in kwargs.keys():
            self.center = kwargs['center']
        else:
            self.center = [0,0,0]
        
    def inside(self, position):
        ''' This checks if a given point is inside of the shape
        '''
        x = position[0]
        y = position[1]
        z = position[2]
        inside = True
        if ((z < (self.center[2] - self.h / 2))
        | (z > (self.center[2] + self.h / 2))
        | (np.sqrt(x**2  + y**2) > self.r**2)):
            inside = False
        else:
            inside = True
        return inside
    
    def getRandPosition(self):
        ''' This returns a random position inside the disc
        '''
        theta = rand()*2*np.pi - np.pi
        r = rand() * self.r
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = rand()*self.h - self.h / 2 + self.center[2]
        return np.array([x,y,z])
    
    def getSlicePosition(self, depth, height):
        ''' This function generates a random position inside of a 'slice'
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
    
    def findIntersect(self, position, direction):
        """This finds the intersection point between a line (characterized by
        a position and a direction), and a disc.
        Parameters
        ----------
        position : 3-by-1 ARRAY of FLOATS
            The position vector relative to the center of the boundary.
        direction : 3-by-1 ARRAY of FLOATS.
            The velocity vector in the coordinate system of the boundary.

        Returns
        -------
        intersect : 3-by-1 ARRAY fo FLOATS
            The coordinates where the electron intersects the Disc.
        """

        v = direction / np.linalg.norm(direction)
        
        c = [0,0,1]

        dot = np.dot(v,c)
        
        angle = np.arccos(dot)
        
        if angle == 0: # This checks if the electron is parallel to disc.
            new_z = self.center[2] + self.h /2
            new_p = position
            new_p[2] = new_z
        else:
            p = position
            ''' First find where the line intersects the circle, i.e. the 
            projection of the cylinder along its axis.
            '''
            t1 = (1/(v[1]**2 + v[0]**2) 
                * (-v[1]*p[1] - v[0]*p[0] 
                    + np.sqrt(self.r**2*(v[1]**2+v[0]**2) 
                        - (v[1]*p[0])**2 + (2*v[1]*v[0]*p[1]*p[0])
                        - (v[0]*p[1])**2)))
            t2 = (-1/(v[1]**2 + v[0]**2) 
                * (v[1]*p[1] + v[0]*p[0] 
                    + np.sqrt(self.r**2*(v[1]**2+v[0]**2) 
                        - (v[1]*p[0])**2 + (2*v[1]*v[0]*p[1]*p[0])
                        - (v[0]*p[1])**2)))
            ''' Pick the t that is positive, because we want only the intersection
            in the direction the particle is traveling
            '''
            t = max(t1,t2) 
            new_p = p + t*v
            ''' Then check if the particle is intersecting the caps of the disc
            '''
            if (new_p[2] > (self.center[2] + self.h / 2)):
                t = ((self.center[2]+self.h/2)-p[2]) / v[2]
                new_p = p + t*v
            elif (new_p[2] < (self.center[2] - self.h / 2)):
                t = ((self.center[2]-self.h/2) - p[2]) / v[2]
                new_p = p + t*v        
        intersect = new_p
        
        return intersect

#%%
        
if __name__ == '__main__':
    # Instantiate a Disc of radius 5000 and height 100
    d = Disc(5000,200)  
    d_xyz = []
    for i in range(3500):
        p = d.findIntersect(d.getRandPosition(),d.getRandPosition())
        d_xyz += [p]
    d_xyz = np.array(d_xyz)
    
    
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111, projection='3d')
    xs = d_xyz[:,0]
    ys = d_xyz[:,1]
    zs = d_xyz[:,2]
    ax.scatter(xs, ys, zs)
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    
    #%%
    
    for i in range(100):
        d_xyz += [d.getSlicePosition(0.0001,0.01)]
    d_xyz = np.array(d_xyz) 
    d_xyz1 = []
    for i in range(100):
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
    
     
    s = Sphere(7000)
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
    
    fig = plt.figure(figsize=(10,10))
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
    