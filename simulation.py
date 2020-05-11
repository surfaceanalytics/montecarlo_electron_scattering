# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:49:02 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
import scipy.stats
from numpy.random import standard_cauchy as cauchy
from shapes import Sphere, Disc
from poisson_distance import Poisson
from generate_angle import AngleDist
from rotation import length, rotate
from scatterer import Scatterer
from electron import Electron
import pickle

#%%

class Simulation():
    def __init__(self, r):
        ''' The simulation class is for running a Monte Carlo simulation of 
        many electrons travelling through a scattering medium. The simulation
        is run inside a 3D volume, called the Boundary (it is a sphere). 
        The scattering medium, and the source of the electrons is a volume
        called the Source (it is a cylinder at the center of the Boundary).
        Electrons are generated at random positions within the Source. They 
        have an initial kinetic energy, and are given a random initial 
        direction. A distance is randomly chosen from an exponential 
        distribution, which is defined by the Total Cross Section (total_xsect).
        Electrons can experience two kinds of scattering event: 1. elastic 
        and 2. inelastic scattering. Both elastic and inelastic processes have
        their own corss sections. The total cross section is the sum of elastic
        and inelastic. After choosing a random distance from the exponential
        distribtuion, it is randomly chosen whether the collision is elastic
        or inelastic (this is based on the relative cross sections of the two
        processes). Then, depending on the scatteing type, a random angle is
        chosen from an angular spread distribution. Based on the chosen angle,
        a new velocity is calculated (note: if inelastic scttering, the kinetic
        energy decreases each event). 
        The electron's new position is determined from the random distance 
        divided by the previous velocity. If the new position is inside the 
        volume of the source, then the move is accepted. The position and 
        velocity are updated. A counter for number of times elastically
        or inelastically scattered is incremented.
        If the electron's new position is outside of the Source, then the move 
        is rejected, and the electron is projected onto the Boundary by
        extending its prevous velocity to the Boundary.
        When an electron hts the boundary, its simulation is over, and the
        simultation moves on to calculate the next electron. If the electron
        does not intersect the Boundary ofter a user defined number of
        inelastic scattering events, then the simulation of the electon is over
        and the simulation moves to a new electron.
        This process is repeated for as many electrons as desired.
        '''
        self.boundary = Sphere(r)
        self.intersections = []
        self.intersected = False
        self.acceptance_angle = 30 / 180 * np.pi
        self.collected = []
        self.depth = 0 # these parameters are used if one wants to generate 
        # electrons piece-wise inside the Disc. Then the depth inside the
        # Disc where the electrons are generated is given by self.depth
        # and the thickness of the piece is given by self.height
        self.height = 0.01
        self.n_event_cutoff = 20 # max number of inel events per electron
        self.initial_KE = 1000 # initial kinetic energy of electrons
    def addSource(self,r,h):
        ''' This function generates a source of electrons. It has the shape of
        a Disc, with radius r and heigh h given in nanometers'''
        self.source = Disc(r,h)
        
    def addScatteringMedium(self, shape):
        self.scatteringMedium = shape
        
    def addScatterer(self, density, inel_factor, inel_exp, el_factor,el_exp, Z):
        ''' The scatterer represents a medium containing atoms or molecules
        that can scatter electrons. The Scatterer has a density in atoms/nm^3,
        it has an inelastic and elastic scattering cross section factor.
        The factor is used to approximate the cross section as a function of
        kinetic energy. The factors need to be determined by fitting data from
        NIST to a 1/E^x function. The scatterer also has two angular spread 
        functions: one for elastic and one for inelastic scattering. The
        parameters el_angle and inel_angle represend the FWHM of the angular
        spread function (ni radians).
        '''
        self.scatterer = Scatterer(density, inel_factor, inel_exp, 
                                   el_factor,el_exp, Z)
        self.scatterer.setXSect(self.initial_KE)

    def createElectron(self):
        ''' This step is performed at the beginning of each electron simulation
        It generates an initial position and direction inside the Source.
        The electron carries all of its information (i.e. position, velocity,
        time, n_times elastic and inelastically scattered) in an arrary called
        e.vector.
        '''
        self.e = Electron(self.source) # instantiate electron
        self.e.kinetic_energy = self.initial_KE # overwrite the default KE value
        self.e.initCoords(self.depth, self.height) # initialize random position
        # for electron inside the source
        self.e.initVelocity() # initialize a random velocity vector
        self.results = np.array([np.copy(self.e.vector)]) # create place to
        # keep the results
        self.scatterer.setXSect(self.initial_KE) #change cross sections
        # of the scatterer so that they correspond to the electron's KE
        
    def _step(self):
        ''' This function simulations one 'step' in the simulation of an 
        electron's path. It calls a method from the Scatterer, called Scatter()
        which determines the kind of scattering event (elastic or inelastic), 
        the distance the electron travels before the scattering event, the new
        direction of the electron relative to its old direction, in polar coord
        (radians)
        '''
        # First get the distance electron travels before being scattered, as
        # well as the kind of scattering event and the angle of scattering
        kind, d, theta, phi = self.scatterer.Scatter()
        # The new time is the distance travelled, divided by the length of
        # the velocity vector.
        velocity = self.e.vector[3:6]
        position = self.e.vector[0:3]
        time = d / np.linalg.norm(velocity)

        # Determine the new position
        new_position = (position + (velocity * time)) 
        
        if self.scatteringMedium.inside(new_position):
            # Check if next potential scattering event is inside or outside
            # source. If it is still inside the scattering medium, then 
            # the electron is scattered once again
            self.e.vector[0:3] = new_position
            # Determine the new velocity after scattering
            # The velocity vector needs to be transformed into a different
            # coordinate system to perform the rotation.
            new_velocity = rotate(velocity, theta, phi)
            self.e.vector[3:6] = new_velocity
            # Check the kind of scattering (elastic or inelastic) and increment 
            # the count accordingly
            if kind == 'elastic':
                # increment the elastic scattering count
                self.e.vector[8]+=1
            else:
                # increment inelastic scattering count
                self.e.vector[7]+=1
                self.e.changeSpeed(self.scatterer.getDeltaKE())
                self.scatterer.setXSect(self.e.kinetic_energy)
            # Determine the new time
            self.e.vector[6] += time
            self.e.pathlength += d
            
        else: 
            # This condition is run, if the electron's new position 
            # is not inside the scattering medium. Then the scattering event 
            # is rejected and the electron is projected onto the Sphere.
            
            # To determine path length within scattering medium
            # first find the interesection with the surface of the scatterer
            scat_inter = self.scatteringMedium.findIntersect(position, 
                                                             velocity)
            # Get the distance from intersection to previous position
            d = np.linalg.norm(self.e.vector[:3] - scat_inter)
            
            # Append this distance to the pathlengths
            self.e.pathlength += d
            
            # Then find coordinates where electron will intersect the boundary
            intersect = self.boundary.findIntersect(position, velocity)
            # Set the electron's current position coordinates to intersection
            # coordinates
            self.e.vector[0:3] = intersect
            # Get new time
            distance = np.linalg.norm(intersect - position)
            new_time = distance / np.linalg.norm(velocity)
            self.e.vector[6] = new_time
            # Collect all intersected electrons into a list
            self.intersections += [list(self.e.vector)]
            self.intersected = True
        
    def Simulate(self):
        ''' This function runs a complete set of simulation steps on an 
        electron. It starts by generating the electron, then repeatedly
        performs simulation steps until the electron has either intersected
        the Boundary, or it has been scattered the set number of times 
        (n_event_cutoff). The calcuation results are saved in the attribute 
        self.results
        '''
        self.createElectron()
        self.intersected = False # this consdition is changed inside of the 
        # function self._step()
        while ((not self.intersected) # while elec. has not reached boundary
            & (self.e.vector[7] < self.n_event_cutoff)
            & (self.e.kinetic_energy > 0)): # and not scattered too many times
            self._step()
            self.results = np.append(self.results, 
                                     np.array([np.copy(self.e.vector[:])]), 
                                     axis=0)
                 
    def simulateMany(self, n, *args):
        ''' This function will run the simulation for many (n) electrons.
        The argument 'keep all' tells the programs to store all the calculated
        paths of every electron. The result is held in the attrib. 'all_paths',
        and is a list of numpy arrays. Each array has a different shape [x,9].
        If the argument 'start finish' is provided, then the starting propertes
        and final properties of each electron are returned in the attribute
        called 'start_finish'.
        Otherwise, only the final step in the calculation is saved, in the 
        attribute 'intersections'.
        '''
        if 'keep all' in args:
            # this keeps all the intermediate positions and velocities of every
            # electron
            self.all_paths = []
            for step in range(n):
                self.Simulate()
                self.all_paths += [np.copy(self.results)]
        elif 'start finish' in args:
            start = []
            finish = []
            pathlengths = []
            for step in range(n):
                print(step)
                self.Simulate()
                start += [np.copy(self.results[0])]
                finish += [np.copy(self.results[-1])]
                pathlengths += [self.e.pathlength]
            self.start_finish = [np.array(start), np.array(finish), np.array(pathlengths)]
        else:
            for step in range(n):
                self.Simulate() # this keeps the final positions, velocities
                # and scatter count in the attribute called self.results
        self.intersections = np.array(self.intersections)
        
    def densityFromP(self, P, T=300): # provide P in mbar
        R = 138 # in units of nm^3/mbar/K/atom
        density = P / (T * R)
        return density #returned in atoms / nm^3
                
        
#%%
if __name__ == '__main__':
    sim = Simulation(5000) 
    # argument is the radius of the spherical simulation evnrionment
    sim.addSource(6000,6000) # source is a shape that emits electrons
    # arguments are radius and thickness (in nm) of a disc
    sim.addScatterer(59,6,10,47) 
    # arguments are denstiy in atoms/nm^3
    # inel_factor, which determines the inelastic cross section
    # el_factor, which determines the elastic cross section
    # and atomic number
    sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford', 
                                                    energy = sim.initial_KE,
                                                    Z = sim.scatterer.Z,
                                                    param = 0.1), 
                           'inelastic':AngleDist(kind = 'Constant')}
    
    sim.height = 75

#%% 

# very long calculation
if __name__ == '__main__':

    sim.simulateMany(1000, 'start finish')
    
    results = sim.start_finish


#%% 
# plot all intersections with the boundary
    xyz = intersections[:,0:3]
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111, projection='3d')
    xs = xyz[:,0]
    ys = xyz[:,1]
    zs = xyz[:,2]
    ax.scatter(xs, ys, zs)
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show()

    # plot of inelsatic scattering events    
    plt.hist(coll[:,7], bins=20, range=(1,20))
    # plot of elsatic scattering events    
    plt.hist(coll[:,8], bins=20)

#%% 
# plot all intersections with the boundary
    data = coll        
    xyz = data[:,0:3]
    fig = plt.figure(figsize = (6,3))
    ax = fig.add_subplot(111, projection='3d')
    xs = xyz[:,0]
    ys = xyz[:,1]
    zs = xyz[:,2]
    ax.scatter(xs, ys, zs)
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show()

    # plot of inelsatic scattering events    
    plt.hist(coll[:,7], bins=20)
    # plot of inelsatic scattering events    
    plt.hist(intersections[:,7], bins=20)

#%%
# this gets the angular distribution of the collected electrons
    data = coll
    counts = {}
    old_theta = 0
    for theta in np.arange(0, np.pi/2,np.pi/40):
        print(theta)
        counts[round(theta / np.pi * 180,2)]=np.array(
                [i for i in data if ((np.arcsin(i[5] 
                / np.sqrt(i[3]**2 + i[4]**2 + i[5]**2))
                < theta)
                & 
                (np.arcsin(i[5] 
                / np.sqrt(i[3]**2 + i[4]**2 +i[5]**2))
                > old_theta))])
        old_theta = theta
    
    all = {}
    for i in counts:
        u = [j for j in counts[i]]
        all[i] = u
    
    all = np.array([[k,len(v)] for k,v in all.items()])
    
    
    
    unscattered = {}
    for i in counts:
        u = [j for j in counts[i] if ((j[7] == 0) & (j[8] == 0))]
        unscattered[i] = u
    
    un = np.array([[k,len(v)] for k,v in unscattered.items()])
        
    scattered_1 = {}
    for i in counts:
        u = [j for j in counts[i] if (j[7] == 1)]
        scattered_1[i] = u
        
    once = np.array([[k,len(v)] for k,v in scattered_1.items()])
    
    scattered_2 = {}
    for i in counts:
        u = [j for j in counts[i] if (j[7] == 2)]
        scattered_2[i] = u
        
    twice = np.array([[k,len(v)] for k,v in scattered_2.items()])
    
    scattered_3 = {}
    for i in counts:
        u = [j for j in counts[i] if (j[7] == 3)]
        scattered_3[i] = u
        
    thrice = np.array([[k,len(v)] for k,v in scattered_3.items()])
    
    #plt.plot(all[:,0],all[:,1])
    plt.plot(un[:,0],un[:,1])
    plt.plot(once[:,0],once[:,1])
    plt.plot(twice[:,0],twice[:,1])
    plt.plot(thrice[:,0],thrice[:,1])
    
    plt.show
    
    #%%
    
    data = {'n_event_cufoff': sim.n_event_cutoff, 'boudary_radius':sim.boundary.r,
            'source_radius':sim.source.r, 'source_height':sim.source.h,
            'slice_thickness':1, 'slice_start':0, 'slice_stop':100,
            'electrons_per_slice':15000, 'scatterer_density':59,
            'inel_factor':6, 'el_factor':10,'inel_angle':10,'el_angle':300,
            'poisson_distance_cutoff': 1E-5, 'inel_angle_dist':'Cauchy',
            'el_angle_dist':'Cauchy','duration':elapsed,
            'avg_loss_per_event':10,'results':sim.intersections}
    
    filename = 'simulation_12'
    outfile = open(filename,'wb')
    pickle.dump(data,outfile)
    outfile.close()
    

