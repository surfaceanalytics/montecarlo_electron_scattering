# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:23:29 2020

@author: Mark
"""
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt

try:
    from .generate_angle import AngleDist

except ImportError as e:
    if str(e) == "attempted relative import with no known parent package":
        from generate_angle import AngleDist
    else:
        raise
#%%

class Shape():
    """A class to represent a shape."""

    def __init__(self):
        """Construct object."""
        pass

    def inside(self, position: np.ndarray) -> bool:
        """Determine if a point is inside the shape or not.

        Parameters:
        ----------
            position: array-like object
                Should have format [x,y,z].
        Returns:
        --------
            inside: bool
                if True, then point is inside the shape. Otherwise it is outside
                the shape.
        """
        pass

    def getIntersection(self, position: np.ndarray, direction:np.ndarray) -> np.ndarray:
        """Find the intersection between a line and the shape.

        The line is defined by a position vector and a direction vector.
        The intersection is a position vector.

        Parameters
        ----------
        position : 3-by-1 ARRAY of FLOATS
            The position vector relative to the center of the boundary.
        direction : 3-by-1 ARRAY of FLOATS.
            The velocity vector in the coordinate system of the boundary.

        Returns
        -------
        intersect : 3-by-1 ARRAY fo FLOATS
            The coordinates where the electron intersects the Shape.
        """
        pass

    def getRandPosition(self):
        """Pick a random position from inside the shape.

        Parameters:
        ----------
            None
        Returns:
        -------
            position: array-like
                Shape [1x3] with the format [x,y,z].
        """
        pass

class Sphere():
    """Sphere is a class that represents the simulation environment.

    Electrons
    are generated form a 'source' shape inside of the sphere, and eventually
    intersect the sphere during the simulation.
    The Sphere class has a method to generate x and y positions, given a z
    It has a method that checks if a point[x,y,z] is inside the sphere. This
    method returns a boolean.
    It also has a method to return the intersection of a vector with the sphere.
    The method takes a vector, that represents direction, as input, and returns
    x,y,z as outputs, representing the intersection with the sphere.
    """

    def __init__(self,r:float):
        """Construct the object."""
        self.r = r
        self.theta = AngleDist(kind = 'Theta')
        self.phi = AngleDist(kind = 'Phi')

    def inside(self, position: np.ndarray) -> bool:
        """Check if a given point is inside of the shape."""
        x = position[0]
        y = position[1]
        z = position[2]
        inside = True
        if np.sqrt(x**2 + y**2 + z**2) > self.r:
            inside = False
        else:
            inside = True
        return inside

    def getIntersection(self, position: np.ndarray, direction: np.ndarray) -> np.ndarray:
        """Find the intersection between a line and the shape.

        The line has a psition and a direction.The function returns
        the coordinates of the intersection.
        """
        v = direction / np.linalg.norm(direction)
        p = position
        dot = np.dot(v,p)
        d1 = -dot + np.sqrt(dot**2-(np.linalg.norm(p)**2-self.r**2))
        d2 = -dot - np.sqrt(dot**2-(np.linalg.norm(p)**2-self.r**2))
        d = np.max([d1,d2])
        intersect = p + (v*d)
        return intersect

    def getRandPosition(self) -> np.ndarray:
        """Return a random position from inside the shape."""
        theta = self.theta.getAngle()
        phi = self.theta.getAngle()
        r = rand() * self.r
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        position = np.array([x,y,z])
        return position

class Disc():
    """The Disc class is used as the shape where electrons are generated.

    It can also be used as the scattering medium. It is defined by the
    attributes r (radius in nm) and h (height in nm).
    """

    def __init__(self, r, h, **kwargs):
        """Construct object.

        Parameters:
        ----------
            r: float
                The radius of the Disc in nm
            h: float
                The height of the Disc in nm
            **kwargs: dict
                center: 1x3 array of floats
                    The Center coordinates of the Disc.
        """
        self.r = r
        self.h = h
        if 'center' in kwargs.keys():
            self.center = kwargs['center']
        else:
            self.center = [0,0,0]

    def inside(self, position: np.ndarray) -> bool:
        """Check if a given position is inside of the shape."""
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

    def getRandPosition(self) -> np.ndarray:
        """Return a random position inside the shape."""
        theta = rand()*2*np.pi - np.pi
        r = rand() * self.r
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = rand()*self.h - self.h / 2 + self.center[2]
        return np.array([x,y,z])

    def getSlicePosition(self, depth:float, height:float) -> np.ndarray:
        """Generate a random position inside of a 'slice'.

        The slice represents a part of the Source disc. It has the same radius
        as the Source, but it has its own height and depth below the top
        surface of the disc.
        """
        theta = rand()*2*np.pi - np.pi
        r = rand() * self.r
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = rand()*(-height) - depth
        return np.array([x,y,z])

    def getIntersection(self, position: np.ndarray, direction: np.ndarray) -> np.ndarray:
        """Find the intersection between a line and the shape.

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

# %% Plot randomly generated points inside a disc.
if __name__ == '__main__':

    d = Disc(5000,200)
    d_xyz = []
    for i in range(3500):
        p = d.getIntersection(d.getRandPosition(),d.getRandPosition())
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


# %% Plot randomly generated positions on a sphere.
    s = Sphere(7000)
    xyz_int = []
    for i in range(2000):
        position = [rand()-0.5,rand()-0.5,rand()-0.5]
        direction = [rand()-0.5,rand()-0.5,rand()-0.5]
        xyz_int += [s.getIntersection(position, direction)]
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
