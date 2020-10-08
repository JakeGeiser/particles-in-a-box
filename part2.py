# -*- coding: utf-8 -*-

#Created on Thu Sep 27 10:50:25 2018
#Completed on Mon Oct 22
#@author: Jake Geiser

#PY 525 hw3_1

import numpy as np
import matplotlib.pylab as py
import copy as cp
### Solving for hw3_2 lattice structure ###
### Useful constatns
N = 16 #16 # number of particles
L = 4. # length of the box
m = 1. # mass of the particles
dt = 0.001 # time step or interval between steps
time = 20. # desired run time in seconds
window = 5
nts = int(time/dt) # how many discrete time steps will occur
### Storing arrays
x1 = np.zeros(N)
y1 = np.zeros(N)
x2 = np.zeros(N)
y2 = np.zeros(N)
x3 = np.zeros(N)
y3 = np.zeros(N)
### Set initial conditions
j = 0.
k = 0.
for i in range(N):
    if k == L:
        k -= L
        j += 1
    x1[i] = k+0.5 ; y1[i] = j+0.5
    k += 1
for it in range(N):
    rand = np.random.randint(0,2)
    rand2 = np.random.choice([-1,1])
    if rand == 0:
        x2[it] = x1[it] + 0.0001* dt * rand2
        y2[it] = y1[it]
    if rand == 1:
        x2[it] = x1[it]
        y2[it] = y1[it] + 0.0001* dt * rand2

### Useful functions
def Force(dx,dy): # dx and dy distance from paricle i - particle i+1
    x = 12/((dx**2+dy**2)**7)*dx - 6/((dx**2+dy**2)**4)*dx
    y = 12/((dx**2+dy**2)**7)*dy - 6/((dx**2+dy**2)**4)*dy
    return [4.*x,4.*y]
def R(xt,yt,xold,yold,Fx,Fy): # Calculate the new positions
    xnew = 2*xt - xold + Fx/m*dt**2.
    ynew = 2*yt - yold + Fy/m*dt**2.
    return [xnew,ynew]
def PBC(px1,py1,px2,py2,length): # Calculate the nearest 'version' of other particles
    rtemp = 2.*length
    rmin = [1,1] # place holder
    for k in [-1.,0.,1.]:
        for l in [-1.,0.,1.]:
            xb = px2 + length * k
            yb = py2 + length * l
            x = px1 - xb
            y = py1 - yb
            if np.sqrt(x**2. + y**2.) < rtemp:
                rmin = [x,y] # error is taking place here and at return?
                rtemp = np.sqrt(x**2. + y**2.)
    return rmin
def chk(px3,py3):
    if px3 < 0: px3 = px3%L + L
    if px3 > L: px3 = px3%L
    if py3 < 0: py3 = py3%L + L
    if py3 > L: py3 = py3%L
    return px3,py3
### Carry out particle motion
py.figure()
for it in range(nts):
    if it%window == 0:
        py.clf() 
        py.plot(x2,y2,'o',markersize=20) # Plot in loop to serve as animation
        py.xlim(0,L)
        py.ylim(0,L)
        py.pause(0.001)
    for t in range(N):
        Fxs = []
        Fys = []
        for be in range(0,t):
            temp = PBC(x2[t],y2[t],x2[be],y2[be],L) # calculate ranges for particle t
            if np.sqrt(temp[0]**2.+temp[1]**2.) < 3.1:
                temp2 = Force(temp[0],temp[1]) # calculate forces for these ranges
                Fxs += [temp2[0]]
                Fys += [temp2[1]]
        for af in range(t+1,N):
            temp = PBC(x2[t],y2[t],x2[af],y2[af],L) # calculate ranges for particle t
            if np.sqrt(temp[0]**2.+temp[1]**2.) < 3.1:
                temp2 = Force(temp[0],temp[1]) # calculate forces for these ranges
                Fxs += [temp2[0]]
                Fys += [temp2[1]]  
        Rtemp = R(x2[t],y2[t],x1[t],y1[t],sum(Fxs),sum(Fys)) # new positions stored
        # Check to see if positions need to be fixed
        x3[t] = Rtemp[0]
        y3[t] = Rtemp[1]
        x3[t], y3[t] = chk(x3[t],y3[t])
    # Shift information
    x1 = cp.deepcopy(x2)
    y1 = cp.deepcopy(y2)
    x2 = cp.deepcopy(x3)
    y2 = cp.deepcopy(y3)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    