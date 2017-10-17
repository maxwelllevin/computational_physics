# -*- coding: utf-8 -*-
"""
Numerical Solutions to Differential Equations

Created on Mon Oct 16 15:02:51 2017

@author: Maxwell Levin
"""

import numpy as np
import pylab as pl

# Example for dX(t)/dt = -x**3 + sint with I.C. X(0) = 1


# ======= Euler's Method ======== #

def f(X,T):
    return - X**3 + np.sin(T)

a = 0                       # t0
b = 10                      # tf
N = 100                     # Number of subdivisions
h = np.abs(b-a)/N           # Step size
x = 1                       # Intial Conditions: X(0) = 1
xPoints = []                # Array to store x values (y-axis on plot)
tPoints = np.arange(a,b,h)  # Array to store t values (x-axis on plot)

for t in tPoints:
    xPoints.append(x)       # Add x value to xPoints
    x = x + h*f(x,t)        # Update x value for our differential equation

pl.plot(tPoints, xPoints, 'b', label='Euler')   # Plot our results


# ======== 4th-Order Runge-Kutta Method ======== #
""" Here we apply the 4th-Order Runge-Kutta Method to approximate
    the solution to a first order ODE given by 
    dX/dt = -X^3 + sin(t)
    with the initial condition that X(t0) = X(0) = 1
    The error of the 4th-Order Runge-Kutta Method is O(h^3)
"""


x = 1                                   # Initial Conditions: X(0) = 1
xPointsRK = []                          # Array to store x values (y-axis on plot)
tPointsRK = np.arange(a,b,h)            # Array to store t values (x-axis on plot)

for t in tPoints:
    k1 = h*f(x,t)                       # Set k1 as defined by 4th Order R-K
    k2 = h*f(x + k1/2, t + h/2)         # Set k2
    k3 = h*f(x + k2/2, t + h/2)         # Set k3
    k4 = h*f(x + k3, t + h)             # Set k4
    xPointsRK.append(x)                 # Add x value to xPoints
    x = x + (k1 + 2*k2 + 2*k3 + k4)/6   # Update x value for our differential equation

pl.plot(tPointsRK, xPointsRK, 'r', label='Runge-Kutta')   # Plot our results
pl.legend()
pl.show()


# ========= Application of R-K Method ========= #

""" We seek the solution to the second-order ODE: 
    ml * d/dt(dTheta/dt) = -mgl * sin(Theta) + al(dTheta/dt) + f(t)
    where m is the mass on the pendulum, l is the length of the 
    pendulum, g is the acceleration due to gravity, a is the damping
    constant, and f(t) is the forcing term.
"""


g = -9.8    # Acceleration due to gravity
d = -1      # Damping constant
m = 1       # Mass of pendulum
l = 1       # Length of pendulum
h = 1e-2    # Set our depth/step size

def forcing(s):   
    return np.sin(s)

def derOmega(w, theta_, t):
    #return theta_
    return -g*np.sin(theta_) + d*w/m + forcing(t)/(m*l)
    
def derTheta(w, theta_, t):
    return w

omega = 1.0                     # Initial Conditions: X(0) = 1
theta = 0.0                     # Set theta initial to be 0
omegaPoints = []                # Array to store omega values (y-axis on plot)
thetaPoints = []                # Array to store theta values 
tPoints = np.arange(a,b,h)      # Array to store time  values (x-axis on plot)


for t in tPoints:
    k1 = h*derOmega(omega, theta, t)                # Set k1 as defined by 4th Order R-K
    k2 = h*derOmega(omega + k1/2, theta, t + h/2)   # Set k2 for omega
    k3 = h*derOmega(omega + k2/2, theta, t + h/2)   # Set k3
    k4 = h*derOmega(omega + k3, theta, t + h)       # Set k4
    
    m1 = h*derTheta(omega, theta, t)                # Set m1 as defined by 4th Order R-K
    m2 = h*derTheta(omega + m1/2, theta, t + h/2)   # Set m2 for omega
    m3 = h*derTheta(omega + m2/2, theta, t + h/2)   # Set m3
    m4 = h*derTheta(omega + m3, theta, t + h)       # Set m4
    
    omegaPoints.append(omega)                       # Add omega value to omegaPoints
    thetaPoints.append(theta)                       # Add theta value to thetaPoints
    
    omega = omega + (k1 + 2*k2 + 2*k3 + k4)/6       # Update omega value
    theta = theta + (m1 + 2*m2 + 2*m3 + m4)/6       # Update theta value
    
    while (theta >= 2*np.pi):                       # Ensure that theta wraps properly
        theta -= 2*np.pi
    while (theta < 0):
        theta += 2*np.pi
   

pl.plot(tPoints, omegaPoints)   # Plot our results
pl.xlim([0,10])                 # Set our domain for the plot
pl.show()                       # Display the plot
