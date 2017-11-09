# -*- coding: utf-8 -*-
"""
Shooting Method 

Created on Wed Nov  8 15:19:49 2017

@author: Maxwell Levin
"""

import numpy as np
import pylab as pl

h= 1e-2


tPts = np.arange(0,10,h)


# Equations of Motion
def dx(x, v, t):
    return v
def dv(x, v, t):
    return -9.8

# Solve the differential equation
def solve(x,v):
    xPts = []
    vPts = []
    for t in tPts:
        k1 = h*dx(x, v, t)
        n1 = h*dv(x, v, t)
        k2 = h*dx(x+0.5*k1,v+0.5*n1,t+0.5*h)
        n2 = h*dv(x+0.5*k1,v+0.5*n1,t+0.5*h)
        k3 = h*dx(x+0.5*k2,v+0.5*n2,t+0.5*h)
        n3 = h*dv(x+0.5*k2,v+0.5*n2,t+0.5*h)
        k4 = h*dx(x+k3, v+n3, t+h)
        n4 = h*dv(x+k3, v+n3, t+h)
        x += (k1+2*k2+2*k3+k4)/6
        v += (n1+2*n2+2*n3+n4)/6
        xPts.append(x)
        vPts.append(v)
    return xPts


# Guess Initial Conditions
def binarySearch():
    v1 = -100
    v2 = 100
    error = 1000
    delta = 1e-6
    while error > delta:
        vp = (v1 + v2)/2
        if ( solve(0,vp)[-1] * solve(0,v1)[-1] > 0):
            v1 = vp
        else:
            v2 = vp
        error = abs( solve(0,v2)[-1] - solve(0,v1)[-1])
    print("The initial velocity required is", (v1 + v2)/2, "m/s")


#==================================================================#

binarySearch()

xPts = solve(0,49)

pl.plot(tPts,xPts)
pl.xlabel("Time")
pl.ylabel("X(t)")