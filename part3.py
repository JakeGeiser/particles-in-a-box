# -*- coding: utf-8 -*-

#Created on Thu Sep 27 10:50:25 2018
#Completed on Mon Oct 22
#@author: Jake Geiser

#PY 525 hw3_1

import numpy as np
import matplotlib.pylab as py
import copy as cp
### Useful constatns
N = 16 #16 # number of particles
L = 4. # length of the box
m = 1. # mass of the particles
dt = 0.0001 # time step or interval between steps
time = 2.5 # desired run time in seconds
window = 25
nts = int(time/dt) # how many discrete time steps will occur
### Storing arrays
x3 = np.zeros(N)
y3 = np.zeros(N)
r_t = np.zeros(int(time/dt))
v_t1 = []
### Set initial conditions from previous problem
x1 = np.array([0.12692186,1.14655152,2.12135533,3.13075539,0.71435339,1.7078084,
2.76639917,3.73238345,0.24519842,1.22906841,2.21510929,3.25110702,
0.92330001,1.90003208,2.83814643,3.81833175])
y1 = np.array([0.43337238,0.66226779,0.48126928,0.44449269,1.58228456,1.49729541,
1.48091992,1.45707639,2.52132384,2.47500839,2.40250353,2.44809026,
3.56358399,3.44055439,3.39514417,3.44174923])
"""x2 = np.array([0.12813631,1.14474598,2.12813627,3.13452539,0.7121143,1.70428552,
2.77176138,3.72894095,0.24511683,1.22839879,2.21769505,3.24973686,
0.91858518,1.89690116,2.83374492,3.82396372])
y2 = np.array([0.43685195,0.65948376,0.48362289,0.44397257,1.5840192,1.49571746,
1.48045777,1.44377887,2.51956455,2.46973809,2.4122407,2.44968079,
3.57012440,3.43870932,3.39977178,3.43913385])"""
x2 = cp.deepcopy(x1)
y2 = cp.deepcopy(y1)

### Useful functions
def Force(dx,dy): # dx and dy distance from paricle i - particle i+1
    x = 12/((dx**2+dy**2)**7)*dx - 6/((dx**2+dy**2)**4)*dx
    y = 12/((dx**2+dy**2)**7)*dy - 6/((dx**2+dy**2)**4)*dy
    return [4.*x,4.*y]
def R(xt,yt,xold,yold,Fx,Fy,T_inc): # Calculate the new positions
    x = 2*xt - xold + Fx/m*dt**2.
    y = 2*yt - yold + Fy/m*dt**2.
    xnew = (x-xt)*(T_inc/100 + 1) + xt
    ynew = (y-yt)*(T_inc/100 + 1) + yt
    return [xnew,ynew]
def PBC(px1,py1,px2,py2,length): # Calculate the nearest ’version’ of other particles
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
T = 0.
for it in range(nts):
    if it%window == 0:
        py.clf()
        py.plot(x2,y2,'o',markersize=20) # Plot in loop to serve as animation
        py.xlim(0,L)
        py.ylim(0,L)
        py.pause(0.001)
    if (it+1)%1000 == 0:
        T += 0.05
        print(T, it*dt)
    rexp = []
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
        for RR in range(t+1,N):
            temp = PBC(x2[t],y2[t],x2[RR],y2[RR],L)
            rexp += [abs(temp[0]**2.+temp[1]**2.)]
        rexp += []
        Rtemp = R(x2[t],y2[t],x1[t],y1[t],sum(Fxs),sum(Fys),T) # new positions
        # Check to see if positions need to be fixed
        x3[t] = Rtemp[0]
        y3[t] = Rtemp[1]
        x3[t], y3[t] = chk(x3[t],y3[t])
    r_t[it] = sum(rexp)/(len(rexp))
    # Shift information
    x1 = cp.deepcopy(x2)
    y1 = cp.deepcopy(y2)
    x2 = cp.deepcopy(x3)
    y2 = cp.deepcopy(y3)
    if it*dt > 0.25 and it*dt < 0.4:
        for t in range(N):
            dist1 = np.sqrt(x1[t]**2. + y1[t]**2.)
            dist2 = np.sqrt(x2[t]**2. + y2[t]**2.)
            v_t1 += [(dist2-dist1)/dt]
py.figure()
rang = np.linspace(0,time,nts)
py.plot(rang,r_t,'.')
X = np.linspace(-3,3,250)
py.figure()
y,x, _ = py.hist(v_t1,bins='auto',range=(-3,3),normed=True)
py.plot(X,np.exp(-X**2./2.)*y.max(),'r')
py.figure()
y,x, _ = py.hist(v_t2,bins='auto',range=(-3,3),normed=True)
13
py.plot(X,np.exp(-X**2./2.)*y.max(),'r')
py.show()

    
    

