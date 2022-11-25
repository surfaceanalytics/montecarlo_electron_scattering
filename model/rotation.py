# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 10:58:22 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt

#%%

def length(vector):
    """Get the length of a vector."""
    v = vector
    l = np.sqrt(np.sum([i **2 for i in v]))
    return l

def rotate(vector: np.ndarray, theta: float, phi:float) -> np.ndarray: # Angles should be provided in radians
    """Rotate a vector using polar coordinates.

    This function takes a velocity vector in some coorindate system. Then
    it transforms the coordinate system to a system, where one of the basis
    vectors is parallel to the velocity vector. Then it rotates the velocity
    vector in in the xz plane about the y axis. Then it transforms back to the
    original coordinate system.

    Parameters:
    ----------
        vector: np.ndarray
            The vector that will be rotated.
            Format is [x,y,z]
        theta: float
            The polar angle of rtation in radians
        phi: float
            The azimuthal angle in radians.

    Returns:
        v1: np.ndarray
            The rotated vector.
            Format is [x,y,z]
    """
    """This is the velocity vector."""
    v0 = vector

    """This generates a random vector that
    should not be parallel to v0 (but no guarantees)
    """
    temp = [rand(),rand(),rand()]

    u1 = v0 / length(v0)
    u2 = np.cross(u1,temp)
    u2 = u2 / length(u2)
    u3 = np.cross(v0,u2)
    u3 = u3 / length(u3)
    C = np.stack([u3,u2,u1],axis=1)
    C_inv = np.linalg.inv(C)
    R0 = np.stack([np.array([np.cos(theta),0,-np.sin(theta)]),
              np.array([0,1,0]),
              np.array([np.sin(theta),0,np.cos(theta)])])
    R1 = np.stack([np.array([np.cos(phi),-np.sin(phi),0]),
                   np.array([np.sin(phi),np.cos(phi),0]),
                   np.array([0,0,1])])

    """This is the rotated vector in the original coordinate system."""
    v1 = np.dot(np.dot(R1,np.dot(R0, np.dot(v0,C))), C_inv)
    return v1

#%%
if __name__ == '__main__':

    v0 = np.array([0,0,1])

    theta = 30 /180 * np.pi
    phi = 10 / 180 * np.pi

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    for phi in range(0,90,1):
        v1 = rotate(v0,theta,phi)
        for i in [v0,v1]:
            xs = [0,i[0]]
            ys = [0,i[1]]
            zs = [0,i[2]]
            ax.plot(xs, ys, zs)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.view_init(0,90)

    plt.show()


