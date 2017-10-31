# -*- coding: utf-8 -*-
"""
Solar System Problem

Created on Mon Oct 30 15:20:25 2017

@author: Maxwell Levin
"""
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:26:32 2017

@author: Maxwell Levin
"""


import numpy as np
import pylab as pl


""" We seek the solution to the second-order ODE: 
    d^2r/dt^2 = -GMsun/r^3 in the outward radial direction.
    We split this into x,y,z coordinates and ignore the z-component
    due to cancellation to achieve the system of ODE's:
    dx/dt = Vx,         dVx/dt = -xGMsun/(x^2 + y^2)^3/2,
    dy/dt = Vy,         dVy,dt = -yGMsun/(x^2 + y^2)^3/2,
    where we set -G*Msun = 4pi^2 to model perfectly circular orbit
"""

ts = 0      # Time start
tf = 1      # Time final ( Years )
h = 1e-2    # Step size ( keep (tf-ts)/h within a reasonable size )
alpha = 0   # Set to 0 for perfectly circular orbit  

tPoints = np.arange(ts,tf,h)    # Array to store our times


# System of ODE's:
def F1(x,y,vx,vy):
    return vx
def F2(x,y,vx,vy):
    return vy
def F3(x,y,vx,vy):
    return -4*np.pi*np.pi*x / ((x*x + y*y)**1.5) * (1 + alpha / (vx*vx + vy*vy))
def F4(x,y,vx,vy):
    return -4*np.pi*np.pi*y / ((x*x + y*y)**1.5) * (1 + alpha / (vx*vx + vy*vy))


# Solve our differential equation using fixed initial conditions
# Using the 4th-order runge-kutta method
def earth_runge_kutta():
    # Arrays to store our values
    xPoints = []
    yPoints = []
    vxPoints = []
    vyPoints = []

    # Initial conditions specific to Earth's Orbit
    x = 1
    y = 0
    vx = 0
    vy = 2*np.pi
    
    for t in tPoints:
        # Runge-Kutta constants
        K1 = h*F1(x, y, vx, vy)
        L1 = h*F2(x, y, vx, vy)
        M1 = h*F3(x, y, vx, vy)
        N1 = h*F4(x, y, vx, vy)
        
        K2 = h*F1(x + K1/2, y + L1/2, vx + M1/2, vy + N1/2)
        L2 = h*F2(x + K1/2, y + L1/2, vx + M1/2, vy + N1/2)
        M2 = h*F3(x + K1/2, y + L1/2, vx + M1/2, vy + N1/2)
        N2 = h*F4(x + K1/2, y + L1/2, vx + M1/2, vy + N1/2)        
        
        K3 = h*F1(x + K2/2, y + L2/2, vx + M2/2, vy + N2/2)
        L3 = h*F2(x + K2/2, y + L2/2, vx + M2/2, vy + N2/2)
        M3 = h*F3(x + K2/2, y + L2/2, vx + M2/2, vy + N2/2)
        N3 = h*F4(x + K2/2, y + L2/2, vx + M2/2, vy + N2/2) 
        
        K4 = h*F1(x + K3, y + L3, vx + M3, vy + N3)
        L4 = h*F2(x + K3, y + L3, vx + M3, vy + N3)
        M4 = h*F3(x + K3, y + L3, vx + M3, vy + N3)
        N4 = h*F4(x + K3, y + L3, vx + M3, vy + N3)
        
        # Store values 
        xPoints.append(x)
        yPoints.append(y)
        vxPoints.append(vx)
        vyPoints.append(vy)
        
        # Update values
        x += (K1 + 2*K2 + 2*K3 + K4)/6
        y += (L1 + 2*L2 + 2*L3 + L4)/6
        vx += (M1 + 2*M2 + 2*M3 + M4)/6
        vy += (N1 + 2*N2 + 2*N3 + N4)/6
    
    # Return the points, capture the last values 
    return xPoints + [x], yPoints + [y]
    

# ==================================================================================== #


xPlot, yPlot = earth_runge_kutta()
pl.plot(xPlot, yPlot)
pl.axes().set_aspect('equal', 'datalim')

