# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 11:57:59 2020

@author: Mark
"""

import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
from poisson_distance import Poisson
import pickle


#%%
data = {'xs1':[],'xs2':[],'xs3':[],'effxs1':[],'effxs2':[],'effxs3':[], 'area1':[],'area2':[]}

#%%
def smooth(x,window_len=11,window='hanning'):

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y   

#%%
''' Image you have two scattering types whose probability
desnity functiuons are given by p1 and p2.
The processes have cross sections xs1 and xs2, respectively
Then the total cross section would be xs3.
We generate Poisson PDF's for each of the three cross sections
and we plot them.
'''
xs1 = 3
xs2 = 1
xs3 = xs1+xs2

p1 = Poisson(xs1,1)
p2 = Poisson(xs2,1)
p3 = Poisson(xs3,1)

trials1 = []
trials2 = []
trials3 = []
for i in range(500000):
    x,a = p1.getDistance()
    trials1 += [[x,a]]
    x,a = p2.getDistance()
    trials2 += [[x,a]]
    x,a = p3.getDistance()
    trials3 += [[x,a]]
    
d1 = np.array(trials1)[:,0]
d2 = np.array(trials2)[:,0]
d3 = np.array(trials3)[:,0]

dist = 'Distance'

fig = plt.figure()
ax = fig.add_subplot(1,1,1,)
ax.hist(d3, bins = 100, color = 'orange', range = (0,10))
ax.hist(d2, bins = 100, color = 'red',range = (0,10))
ax.hist(d1, bins = 100, color = 'blue', range = (0,10))

plt.xlabel(dist)
plt.ylabel('Counts')
plt.show()


#%%
# compare the elastic and inelastic distances
''' Now we compare the distances of the type a and type 2 scattering 
distributions. The shorter distance is chosen.
'''
D = np.array([d1,d2,d3])


inel = D.T[np.where(D[0]<D[1], True, False)][:,0]
el = D.T[np.where(D[0]>D[1], True, False)][:,1]
total = np.concatenate([inel, el])

print(len(inel)/(len(el)+len(inel)))
print(len(el)/(len(el)+len(inel)))


dist = 'Distance'

fig = plt.figure()
ax = fig.add_subplot(1,1,1,)
ax.hist(d3, bins = 100, color = 'orange', range = (0,4))
ax.hist(el, bins = 100, color = 'red', range = (0,4))
ax.hist(inel, bins = 100, color = 'blue', range = (0,4))
ax.hist(total, bins = 100, color = 'green', range = (0,4))

plt.xlabel(dist)
plt.ylabel('Counts')
plt.show()

pdf_inel = np.histogram(inel,1000)
pdf_el = np.histogram(el,1000)
pdf3 = np.histogram(d3,1000)
pdf_total = np.histogram(total,1000)
plt.plot(pdf_inel[1][:-1],pdf_inel[0])
plt.plot(pdf_el[1][:-1],pdf_el[0])
plt.plot(pdf3[1][:-1],pdf3[0])
plt.plot(pdf_total[1][:-1],pdf_total[0])

plt.show()

max_inel = pdf_inel[1][np.argmax(smooth(pdf_inel[0],window_len=21,window='flat'))]
max_el = pdf_el[1][np.argmax(smooth(pdf_el[0],window_len=21,window='flat'))]
max_3 = pdf3[1][np.argmax(smooth(pdf3[0],window_len=21,window='flat'))]
max_total = pdf_total[1][np.argmax(smooth(pdf_total[0],window_len=21,window='flat'))]


#%%
data['xs1'] += [xs1]
data['xs2'] += [xs2]
data['xs3'] += [xs3]
data['effxs1'] += [1/max_inel]
data['effxs2'] += [1/max_el]
data['effxs3'] += [1/max_total]
data['area1'] += [len(inel)/(len(el)+len(inel))]
data['area2'] += [len(el)/(len(el)+len(inel))]


plt.scatter(data['xs3'],data['effxs3'])

expected = [data['xs1'][i] / (data['xs1'][i]+data['xs2'][i]) for i in range(len(data['xs1']))]

plt.scatter(expected, data['area1'])

#%%

filename = 'poisson_1'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close()


