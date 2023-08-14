# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 10:56:56 2022

@author: pielsticker
"""
import os
import numpy as np
import matplotlib.pyplot as plt
#%% Show that equation can be converted.
theta = np.arange(0,3.2,0.01)

fig, ax = plt.subplots(figsize=(5, 4), dpi=300)
#beta = [0.1,0.2,2,5]
fontdict = {"size": 15}

def _beta(energy, Z, param):
    beta = param * 5.43 * (Z**(2/3))/energy
    return beta

Z = 47
E = 400
p = 1#0.9

beta = _beta(E,Z,p)
R = ((np.cos(theta)-1)*(1+beta))/(np.cos(theta)-2*beta-1)
h = ax.plot(theta, R)

ax.set_xlim(left=np.min(theta), right=np.max(theta))
ax.set_ylim(bottom=0)
ax.set_xlabel("\u03B8 (radians)",fontdict=fontdict)
ax.set_ylabel("R",fontdict=fontdict)
ax.tick_params(axis="x", labelsize=fontdict["size"])
ax.tick_params(axis="y", labelsize=fontdict["size"])

#ax.set_title()
fig.tight_layout()

fig.savefig("rutherford_invertible.png", bbox_inches="tight")
#%% Sample random values of theta using R in (0,1].
R = np.arange(0,1,0.01)

fig, ax = plt.subplots(figsize=(5, 4), dpi=300)
beta = [0.1,0.2,2,5]
fontdict = {"size": 15}

beta = _beta(E,Z,p)
theta = np.arccos((-1-beta+R+2*R*beta)/(-1-beta+R))
h = ax.plot(R, theta)

ax.set_xlim(left=np.min(R), right=np.max(R))
ax.set_ylim(bottom=0)
ax.set_xlabel("R",fontdict=fontdict)
ax.set_ylabel("\u03B8 (radians)",fontdict=fontdict)
ax.tick_params(axis="x", labelsize=fontdict["size"])
ax.tick_params(axis="y", labelsize=fontdict["size"])

#ax.set_title()
fig.tight_layout()

fig.savefig("rutherford_inverted.png", bbox_inches="tight")