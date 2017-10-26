# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:13:48 2017

@author: Maxwell Levin
"""


""" We seek the solution to the second-order ODE: 
    ml * d/dt(dTheta/dt) = -mgl * sin(Theta) + al(dTheta/dt) + f(t)
    where m is the mass on the pendulum, l is the length of the 
    pendulum, g is the acceleration due to gravity, a is the damping
    constant, and f(t) is the forcing term.
"""

import numpy as np
import pylab as pl

# Example for dX(t)/dt = -x**3 + sint with I.C. X(0) = 1



a = 0                       # t0
b = 100                      # tf
g = 9.8    # Acceleration due to gravity
q = 1/2     # Damping constant
m = 1       # Mass of pendulum
L = 9.8       # Length of pendulum
h = 1e-2    # Set our depth/step size
w_d = 2/3   # Set the frequency of our forcing
f_d = 0.2
h = 1e-3

def forcing(s):   
    return np.sin(s)

def F2(theta_, w, t):
    #return theta_
    return -g*np.sin(theta_)/L -q*w + f_d*forcing(w_d*t) 
    
def F1(theta_, w, t):
    return w



omega = 0.0                     # Initial Conditions: X(0) = 1
theta = 0.2                     # Set theta initial to be 0
omegaPoints = []                # Array to store omega values (y-axis on plot)
thetaPoints = []                # Array to store theta values 
tPoints = np.arange(a,b,h)      # Array to store time  values (x-axis on plot)


for t in tPoints:
    K1 = h*F1(theta, omega, t)
    L1 = h*F2(theta, omega, t)
    K2 = h*F1(theta + K1/2, omega + L1/2, t + h/2)
    L2 = h*F2(theta + K1/2, omega + L1/2, t + h/2)
    K3 = h*F1(theta + K2/2, omega + L2/2, t + h/2)
    L3 = h*F2(theta + K2/2, omega + L2/2, t + h/2)
    K4 = h*F1(theta + K3, omega + L3, t + h)
    L4 = h*F2(theta + K3, omega + L3, t + h)
    
    omegaPoints.append(omega)                       # Add omega value to omegaPoints
    thetaPoints.append(theta)                       # Add theta value to thetaPoints
    
    theta += (K1 + 2*K2 + 2*K3 + K4)/6
    omega += (L1 + 2*L2 + 2*L3 + L4)/6
    
    if (theta < -np.pi):                       # Ensure that theta wraps properly
        theta += 2*np.pi
    if (theta > np.pi):
        theta -= 2*np.pi
   

pl.plot(tPoints, thetaPoints)       # Plot our results
pl.show()
pl.plot(thetaPoints, omegaPoints)   # Plot our results
pl.show()                           # Display the plot