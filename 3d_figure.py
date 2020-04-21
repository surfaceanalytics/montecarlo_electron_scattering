# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 17:39:15 2020

@author: Mark
"""
from environment import Simulation
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
import time
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import cm
#%%

sim = Simulation(5000) # argument is the radius of the spherical simulation evnrionment
sim.addSource(2000,5000) # source is a shape that emits electrons
# arguments are radius and thickness of a disc
sim.addScatterer(59,6,10,1/180*np.pi,0/180*np.pi) # arguments are denstiy in atoms/nm^3

sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Cauchy', width = 1/180*np.pi), 
                           'inelastic':AngleDist(kind = 'Cauchy', width = 50/180*np.pi)}
# inel_factor, which determines the inelastic cross section
# el_factor
# inelastic angular spread in degrees, and elastic angular spread in degrees
sim.depth = 0
sim.height = 75

#%%
# this is to check that the initial velocity directions are spherically distributed

V = []
for i in range(5000):
    sim.createElectron()
    V += [sim.e.vector]
V = np.array(V)

V[:,:3] = V[:,:3]*0

fig = plt.figure(figsize=(6,6))
ax = fig.gca(projection='3d')

# Make the grid
x, y, z = V[:,0], V[:,1], V[:,2]

# Make the direction data for the arrows
u, v, w = V[:,3], V[:,4], V[:,5]

ax.quiver(x, y, z, u, v, w, length=0.05, normalize=True)

plt.show()

lengths = [length(i[3:6]) for i in V]

#%%
#This calculated the trajectories of j electrons, and keeps the entire path
manysteps = []
n_electrons = 500000
n_steps = 40
for j in range(n_electrons):
    print(j)
    sim.createElectron()
    steps = [np.copy(sim.e.vector)]
    for i in range(n_steps):
        sim._step()
        steps +=[np.copy(sim.e.vector)]
    steps = np.array(steps)
    manysteps += [np.copy(steps)]
manysteps = np.array(manysteps)
#%%
#Set color map so that colors are mapped to number of times inelastically scattered

max_z = np.max(manysteps[:,:,7])
start = 0
stop = 1
n_lines = int(max_z)+1
cm_subsection = np.linspace(start, stop, n_lines)
colors = [ cm.plasma(x) for x in cm_subsection]

#%%
# This plots all the paths of the electrons
fig = plt.figure(figsize=(10,10))
xlim = 200
ylim = xlim
zlim = 200
ax = fig.gca(projection='3d')
xLabel = ax.set_xlabel('x distance', linespacing=3)
yLabel = ax.set_ylabel('y distance', linespacing=3)
zLabel = ax.set_zlabel('z distance', linespacing=3)
c = np.array([s[-1,7] for s in manysteps])
for k in manysteps[:10000]:
    V = k
    n_scatter = int(np.max(V[:,7])) # get number of times scattered
    color = colors[n_scatter]
    x, y, z = V[:,0], V[:,1], V[:,2]
    
    # Make the direction data for the arrows
    #u, v, w = V[:,3], V[:,4], V[:,5]
    
    #ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)
    ax.plot3D(x, y, z, c = color)

    ax.set_xlim3d(-xlim, xlim)
    ax.set_ylim3d(-ylim,ylim)
    ax.set_zlim3d(-zlim,zlim)

ax.set_xticks([-xlim, 0, xlim])
ax.set_yticks([-ylim, 0, ylim])
ax.set_zticks([-zlim, 0, zlim])

ax.view_init(0, 90)    
plt.show()



#%%
# This plots all final positions of the electrons, where the color is number of time inel scattered

from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')

fig = plt.figure(figsize=(10,10))
xlim = 5000
ylim = xlim
zlim = xlim
ax = fig.gca(projection='3d')
xLabel = ax.set_xlabel('x distance', linespacing=3)
yLabel = ax.set_ylabel('y distance', linespacing=3)
zLabel = ax.set_zlabel('z distance', linespacing=3)
c = np.array([s[-1,7] for s in manysteps])
for k in manysteps:
    V = k
    n_scatter = int(np.max(V[:,7])) # get number of times scattered
    color = colors[n_scatter]
    # Make the grid
    x, y, z = V[:,0], V[:,1], V[:,2]
    
    # Make the direction data for the arrows
    #u, v, w = V[:,3], V[:,4], V[:,5]
    
    #ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)
    #ax.plot3D(x, y, z, c = color)
    ax.scatter3D(x[-1], y[-1], z[-1], c = color)

    ax.set_xlim3d(-xlim, xlim)
    ax.set_ylim3d(-ylim,ylim)
    ax.set_zlim3d(-zlim,zlim)

ax.set_xticks([-xlim, 0, xlim])
ax.set_yticks([-ylim, 0, ylim])
ax.set_zticks([-zlim, 0, zlim])

ax.view_init(0, 90)    
plt.show()


#%%

positive = 0
above = []
negative = 0
below = []
for i in manysteps:
    if i[-1,2] > 0:
        positive += 1
        above += [i[-1]]
    else:
        negative += 1
        below += [i[-1]]
        
above = np.array(above)
below = np.array(below)

#%%
# This plots the intesection of the electrons that have been scattered n times
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
last_position = manysteps[:,-1,:]
n_times_scattered = 0
once_only = last_position[np.where(last_position[:,7]==n_times_scattered)]
fig = plt.figure(figsize=(10,10))
xlim = 5000
ylim = xlim
zlim = xlim
ax = fig.gca(projection='3d')
xLabel = ax.set_xlabel('x distance', linespacing=3)
yLabel = ax.set_ylabel('y distance', linespacing=3)
zLabel = ax.set_zlabel('z distance', linespacing=3)
c = np.array([s[-1,7] for s in manysteps])
for k in once_only:
    V = k
    x, y, z = V[0], V[1], V[2]
    
    ax.scatter3D(x, y, z, c = colors[n_times_scattered])

    ax.set_xlim3d(-xlim, xlim)
    ax.set_ylim3d(-ylim,ylim)
    ax.set_zlim3d(-zlim,zlim)

ax.set_xticks([-xlim, 0, xlim])
ax.set_yticks([-ylim, 0, ylim])
ax.set_zticks([-zlim, 0, zlim])

ax.view_init(90, 90)    
plt.show()

#%%
# This plots the staring points of electrons that have collided with the sphere
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
# first pick only the electrons that hit the upper half of the sphere
hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]
# then get the number of times these electrons were inelastically scattered
scatter_count= hitting_sphere[:,-1,7]
# then get the initial position of these electrons
first_position = hitting_sphere[:,0,:] 

max_z = np.max(scatter_count)
start = 0
stop = 1
n_lines = int(max_z)+1
cm_subsection = np.linspace(start, stop, n_lines)
colors = [ cm.plasma(x) for x in cm_subsection]

fig = plt.figure(figsize=(10,10))
xlim = 1500
ylim = xlim
zlim = 1500
ax = fig.gca(projection='3d')
xLabel = ax.set_xlabel('x distance', linespacing=3)
yLabel = ax.set_ylabel('y distance', linespacing=3)
zLabel = ax.set_zlabel('z distance', linespacing=3)
c = np.array([s[-1,7] for s in manysteps])
for l, k in enumerate(first_position[:10000]):
    V = k
    x, y, z = V[0], V[1], V[2]
    # color the dots according to the number of times the electron was scattered
    ax.scatter3D(x, y, z, c = colors[int(scatter_count[l])])

    ax.set_xlim3d(-xlim, xlim)
    ax.set_ylim3d(-ylim,ylim)
    ax.set_zlim3d(-zlim,zlim)

ax.set_xticks([-xlim, 0, xlim])
ax.set_yticks([-ylim, 0, ylim])
ax.set_zticks([-zlim, 0, zlim])

ax.view_init(10, 90)    
plt.show()

#%%
# This plots the average number of times being inelastically scattered versus
# depth
# first pick only the electrons that hit the upper half of the sphere
hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]
# then get the number of times these electrons were inelastically scattered
scatter_count = hitting_sphere[:,-1,7]
# then get the initial position of these electrons
first_position = hitting_sphere[:,0,:] 

# change scatter count column 
first_position[:,7] = scatter_count

# find min and max z positions
min_z = np.min(first_position[:,2])
max_z = np.max(first_position[:,2])
step_size = 0.5
steps = np.arange(-50,0,step_size)
avg_depth = []
cts_depth = []

for s in steps:
    avg = np.mean(first_position[np.where((first_position[:,2] < s + step_size)
    &(first_position[:,2] >= s)),7])
    cts = len((first_position[np.where((first_position[:,2] < s + step_size)
    &(first_position[:,2] >= s)),7][0,:]))
    avg_depth += [[s,avg]]
    cts_depth += [[s,cts]]
avg_depth = np.array(avg_depth)
cts_depth = np.array(cts_depth)

# this figure shows the number of times inelastically scattered versus the
# average depth below the surface for electrons being scattered that many times
fig = plt.figure(figsize=(5,5))
ax = fig.gca()
ax.plot(avg_depth[:,0], avg_depth[:,1])
ax.set_xlabel('Depth [nm]', linespacing=3)
ax.set_ylabel('Number of times inelastically scattered', linespacing=3)
plt.show()

# This plot shows the number of electrons reaching the boundary versus depth
fig = plt.figure(figsize=(5,5))
ax = fig.gca()
ax.plot(cts_depth[:,0], cts_depth[:,1])
ax.set_xlabel('Depth [nm]', linespacing=3)
ax.set_ylabel('Number of electrons reaching boundary', linespacing=3)
plt.show()



#%%
# this shows the number of electrons reachng the boundary versus depth for 
# a given number of scattering events
# the plots should have a Poisson distribution

# first pick only the electrons that hit the upper half of the sphere
hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]
# then get the number of times these electrons were inelastically scattered
scatter_count = hitting_sphere[:,-1,7]
# then get the initial position of these electrons
first_position = hitting_sphere[:,0,:] 

# change scatter count 
first_position[:,7] = scatter_count

# This plots shows the number of electrons reaching the boundary for a given 
# n. The plots should follow the Poisson distribution
fig = plt.figure(figsize=(5,5))
ax = fig.gca()

all_distributions = []
for n in range(7):
    single_n = first_position[np.where(first_position[:,7] == n)]
    
    # find min and max z positions
    min_z = np.min(first_position[:,2])
    max_z = np.max(first_position[:,2])
    step_size = 0.5
    steps = np.arange(-20,0,step_size)
    cts_depth = []
    
    for s in steps:
        cts = len((single_n[np.where((single_n[:,2] < s + step_size)
        &(single_n[:,2] >= s)),7][0,:]))
        cts_depth += [[s,cts]]
    cts_depth = np.array(cts_depth)
    all_distributions +=[cts_depth[:,1]]

    
    ax.plot(cts_depth[:,0], cts_depth[:,1])

ax.set_xlabel('Depth [nm]', linespacing=3)
ax.set_ylabel('Number of electrons reaching boundary', linespacing=3)
plt.show()

all_distributions = np.array(all_distributions).T

sums = np.sum(all_distributions, axis=0)
#%%
# This plots the average path length electrons travel versus number of times they 
# are iinelastically scattered (only electrons that hit the detector)
# first pick only the electrons that hit the upper half of the sphere
hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]

paths = []
path_lengths = []
for i in hitting_sphere:
    path = np.linalg.norm(i[1:,:3]-i[:-1,:3], axis = 1)
    paths += [path[np.where((path != 0) & (path < 4000))]]
    length = np.sum(np.where((path != 0) & (path < 4000)))
    n_event = i[-1,7]
    path_lengths += [[n_event,length]]
path_lengths = np.array(path_lengths)

avg_distances = []
for i in range(int(np.max(path_lengths[:,0]))):
    dist = path_lengths[np.where(path_lengths[:,0]==i)]
    avg_distances += [[i,np.mean(dist[:,1])]]
avg_distances = np.array(avg_distances)

fig = plt.figure(figsize=(8,5))
ax = fig.gca()
ax.scatter(avg_distances[:,0],avg_distances[:,1])
ax.set_xlabel('Number of times scattered', linespacing=3)
ax.set_ylabel('Average path length [nm]', linespacing=3)
plt.show()

#%%
def sphere_intersections(manysteps):
    # first pick only the electrons that hit the upper half of the sphere
    hitting_sphere = manysteps[np.where(manysteps[:,-1,2]>0)]
    # then get the number of times these electrons were inelastically scattered
    scatter_count = hitting_sphere[:,-1,7]
    # then get the final position and velocity of these electrons
    final = hitting_sphere[:,-1,[0,1,2,3,4,5,7]] 
    return final

def get_angle_dist(intersections, start = 0, stop = np.pi/2, step = np.pi/100, **kwargs):
    #intersections = sphere_intersections(manysteps)
    if 'n_scatter' in kwargs.keys():
        n = kwargs['n_scatter']
        intersections = intersections[np.where(intersections[:,-1] == int(n))]
    
    velocities = intersections[:,3:6]
    speeds = np.linalg.norm(velocities,axis=1)
    angles = np.array([np.arctan(velocities[:,1] / velocities[:,0]), np.arccos(velocities[:,2] / speeds)]).T
    angle_dist = []
    r = np.mean(speeds)
    for i in np.arange(start,stop, step):
        count = len(angles[np.where((angles[:,1] > i) & (angles[:,1] <= i + step))])
        area = 2*np.pi*(r**2)*(-np.cos(i + step)+np.cos(i))
        angle_dist += [[i, count / (area/1E+13), count, area]]
    angle_dist = np.array(angle_dist)    
    return angle_dist

angle_dist = get_angle_dist(manysteps, n_scatter = 6)

fig = plt.figure(figsize=(5,5))
ax = fig.gca()
ax.scatter(angle_dist[:,0], angle_dist[:,1])
ax.set_xlabel('Counts / solid angle', linespacing=3)
ax.set_ylabel('Polar Angle [radians]', linespacing=3)
plt.show()

fig = plt.figure(figsize=(5,5))
ax = fig.gca()
ax.scatter(angle_dist[:,0], angle_dist[:,2])
ax.set_xlabel('Counts / angle', linespacing=3)
ax.set_ylabel('Polar Angle [radians]', linespacing=3)
plt.show()

fig = plt.figure(figsize=(5,5))
ax = fig.gca()
for i in range(0,7):
    angle_dist = get_angle_dist(manysteps, n_scatter = i)
    ax.scatter(angle_dist[:,0], angle_dist[:,1])
ax.set_xlabel('Counts / solid angle', linespacing=3)
ax.set_ylabel('Polar Angle [radians]', linespacing=3)
plt.show()