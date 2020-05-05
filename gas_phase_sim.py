# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:49:49 2020

@author: Mark
"""
from simulation import Simulation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from generate_angle import AngleDist
from matplotlib import cm
import pickle
from shapes import Sphere
from analysis import Analysis
#%%

if __name__ == '__main__':
    sim_radius = 800000
    sim = Simulation(sim_radius) 
    # argument is the radius of the spherical simulation evnrionment
    sim.initial_KE = 1100
    sim.addSource(100000,5000) # source is a shape that emits electrons
    # arguments are radius and thickness (in nm) of a disc
    sim.addScatteringMedium(Sphere(sim_radius))
    sim.addScatterer(9.66E-5,0.527,0.85,0.5,1.01,4) 
    # arguments are denstiy in atoms/nm^3
    # inel_factor, and inel_exp which determines the inelastic cross section
    # from a fit to TPP2M and NIST data using sigma = facter / (KE * exp)
    # el_factor, which determines the elastic cross section
    # and atomic number
    sim.scatterer.angle_dist = {'elastic':AngleDist(kind = 'Rutherford', 
                                                    energy = sim.initial_KE,
                                                    Z = sim.scatterer.Z,
                                                    param = 0.06), 
                           'inelastic':AngleDist(kind = 'Constant')}
    #set the high of the source slice
    sim.height = 1

#%% 

# potentially very long calculation
if __name__ == '__main__':

    sim.simulateMany(30000, 'start finish')
    results = sim.start_finish

#%%
if __name__ == '__main__':
    A = Analysis(results)
    A.pathHistogram('show',bins=25, x_limits=(700000,900000), accept_angle =
                    (30,800000))        
    A.areasUnderProfiles('show')

    p = A.profiles
    a = A.areas
    
    l_not_zero = np.histogram(results[2][np.where(results[2] != 0)], bins=100)
    plt.plot(l_not_zero[1][:-1],l_not_zero[0])    
    
#%%
# This plots the average number of times being inelastically scattered versus
# depth
# first pick only the electrons that hit the upper half of the sphere

        
    hitting_sphere = A.results[1]
    acceptance_angle = 30 #degrees
    radius = sim_radius * np.sin(acceptance_angle / 180 * np.pi)
    acceptance = np.array([i for i in hitting_sphere if (((i[0]**2
                  + i[1]**2)
                    < (radius)**2)
                    & (i[2] > 0))])
    
    accepted_lengths = np.array([A.results[2][i] for i in range(len(A.results[2]))
                        if ((A.results[1][i][0]**2
                            + A.results[1][i][1]**2
                            < (radius)**2)
                        & (A.results[1][i][2] > 0))])

    for i in range(len(accepted_lengths)):
        if accepted_lengths[i] == 0:
            end = A.results[1][i][:3]
            start = A.results[0][i][:3]
            length = np.linalg.norm(end-start)
            accepted_lengths[i] = length
        
    path_hist = np.histogram(accepted_lengths, bins=100)
    plt.plot(path_hist[1][:-1], path_hist[0])
    
    fig = plt.figure(figsize=(10,10))
    xlim = 800000
    ylim = xlim
    zlim = xlim
    ax = fig.gca(projection='3d')
    xLabel = ax.set_xlabel('x distance', linespacing=3)
    yLabel = ax.set_ylabel('y distance', linespacing=3)
    zLabel = ax.set_zlabel('z distance', linespacing=3)
    for l, k in enumerate(acceptance[:500]):
        V = k
        x, y, z = V[0], V[1], V[2]

        ax.scatter3D(x, y, z)
    
    for l,k in enumerate(A.results[0][:500]):
        V = k
        x, y, z = V[0], V[1], V[2]

        ax.scatter3D(x, y, z)
        
    ax.set_xlim3d(-xlim, xlim)
    ax.set_ylim3d(-ylim,ylim)
    ax.set_zlim3d(-100000,zlim)
    
    ax.set_xticks([-xlim, 0, xlim])
    ax.set_yticks([-ylim, 0, ylim])
    ax.set_zticks([-zlim, 0, zlim])
    
    ax.view_init(10, 90)    
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
for n in range(12):
    single_n = first_position[np.where(first_position[:,7] == n)]
    
    # find min and max z positions
    min_z = np.min(first_position[:,2])
    max_z = np.max(first_position[:,2])
    step_size = 0.5
    steps = np.arange(-40,0,step_size)
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