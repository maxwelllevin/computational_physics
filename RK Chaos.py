# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:26:32 2017

@author: Maxwell Levin
"""


import numpy as np
import pylab as pl


""" We seek the solution to the second-order ODE: 
    ml * d/dt(dTheta/dt) = -mgL * sin(Theta) + aL(dTheta/dt) + f(t)
    where m is the mass on the pendulum, l is the length of the 
    pendulum, g is the acceleration due to gravity, a is the damping
    constant, and f(t) is the forcing term.
"""

a = 0       # Start of domain
b = 500     # End of domain
h = 2e-1    # Step size ( keep (b-a)/h within a reasonable size )
g = 9.8     # Acceleration due to gravity
q = 1/2     # Damping constant
m = 1       # Mass of pendulum
L = 9.8     # Length of pendulum
h = 1e-2    # Set our depth/step size
w_d = 2/3   # Set the frequency of our forcing
f_d = 1.44  # Strength of the forcing


# forcing function (give it periodic forcing)
def forcing(s):   
    return f_d*np.sin(s)

# Our d[Omega]/dt function
def F2(theta_, w, t):
    return -g*np.sin(theta_)/L - q*w + forcing(w_d*t) 

# Our d[theta]/dt function
def F1(theta_, w, t):
    return w

# Our function that computes the solution to the second order differential
# equation using 4th order Runge-Kutta method
def chaos_runge_kutta(theta, omega):
    omegaPoints = []
    thetaPoints = []
    snapOmega = []
    snapTheta = []
    tPoints = np.arange(a,b,h)
    
    for t in tPoints:
        K1 = h*F1(theta, omega, t)
        L1 = h*F2(theta, omega, t)
        K2 = h*F1(theta + K1/2, omega + L1/2, t + h/2)
        L2 = h*F2(theta + K1/2, omega + L1/2, t + h/2)
        K3 = h*F1(theta + K2/2, omega + L2/2, t + h/2)
        L3 = h*F2(theta + K2/2, omega + L2/2, t + h/2)
        K4 = h*F1(theta + K3, omega + L3, t + h)
        L4 = h*F2(theta + K3, omega + L3, t + h)
        
        omegaPoints.append(omega)                       
        thetaPoints.append(theta)
        
        if (np.abs( w_d*t/(2*np.pi) - np.floor(w_d*t/(2*np.pi))) <= h/10):
            snapOmega.append(omega)
            snapTheta.append(theta)
        
        theta += (K1 + 2*K2 + 2*K3 + K4)/6
        omega += (L1 + 2*L2 + 2*L3 + L4)/6
        
        if (theta < -np.pi):                      
            theta += 2*np.pi
        if (theta > np.pi):
            theta -= 2*np.pi
    
    return thetaPoints, omegaPoints, tPoints, snapTheta, snapOmega
    

# ==================================================================================== #

theta1, omega1, time, snapTheta, snapOmega = chaos_runge_kutta(0.2, 0.0)
theta2, omega2, time, dum1, dum2 = chaos_runge_kutta(0.2+1e-3, 0)

diff = []
for i in range(len(theta1)):
    temp = np.abs(theta1[i] - theta2[i])   
    diff.append(np.log(temp))


# Plot the 'structure' of the differential equation
# by plotting omega vs theta (snapshots are taken under helpful conditions)
pl.plot(snapTheta, snapOmega, "k.")
pl.show()

# Plot theta1 and theta2 on same plot (slightly different initial conditions)
pl.plot(time, theta1, label="theta intial is 0.2")
pl.plot(time, theta2)
pl.title("Theta vs t")
pl.show()

# Plot theta1 vs omega1 and theta2 vs omega2 on same plot
# this helps us see the chaos of the structure
pl.plot(theta1, omega1)
pl.plot(theta2, omega2)
pl.title("Omega vs t")
pl.show()      

# Plot the log of the difference between the theta plots    
pl.plot(time,diff)
pl.show()


