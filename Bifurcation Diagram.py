# -*- coding: utf-8 -*-
"""
Bifurcation Diagram

Created on Fri Oct 27 15:54:42 2017

@author: Maxwell Levin
"""


""" We seek to plot the solution to the differential
    equation we solved in RK Chaos.py, but vary the strength
    of our forcing term linearly so that we can see the
    bifurcation diagram produced by plotting theta against
    the strength of the forcing term. (We expect our 
    diagram to resemble a fractal)
"""

import numpy as np
import pylab as pl



a = 0       # Start of domain
b = 150     # End of domain
h = 1e-1    # Step size ( keep (b-a)/h within a reasonable size )
g = 9.8     # Acceleration due to gravity
q = 1/2     # Damping constant
m = 1       # Mass of pendulum
L = 9.8     # Length of pendulum
h = 1e-2    # Set our depth/step size
w_d = 2/3   # Set the frequency of our forcing


tPoints = np.arange(a,b,h)  # Create the t-axis array

# forcing function (give it periodic forcing)
def forcing(s, f_d):   
    return f_d*np.sin(s)

# Our d[Omega]/dt function
def F2(theta_, w, t, f_d):
    return -g*np.sin(theta_)/L - q*w + forcing(w_d*t, f_d) 

# Our d[theta]/dt function
def F1(theta_, w, t):
    return w

# Our function that computes the solution to the second order differential
# equation using 4th order Runge-Kutta method
def chaos_runge_kutta(theta, omega, f_d):
    omegaPoints = []
    thetaPoints = []
    snapTheta = []
    f_d_list = []
    for t in tPoints:
        K1 = h*F1(theta, omega, t)
        L1 = h*F2(theta, omega, t, f_d)
        K2 = h*F1(theta + K1/2, omega + L1/2, t + h/2)
        L2 = h*F2(theta + K1/2, omega + L1/2, t + h/2, f_d)
        K3 = h*F1(theta + K2/2, omega + L2/2, t + h/2)
        L3 = h*F2(theta + K2/2, omega + L2/2, t + h/2, f_d)
        K4 = h*F1(theta + K3, omega + L3, t + h)
        L4 = h*F2(theta + K3, omega + L3, t + h, f_d)
        omegaPoints.append(omega)                       
        thetaPoints.append(theta)
        if (np.abs( w_d*t/(2*np.pi) - np.floor(w_d*t/(2*np.pi))) <= h/10) and t > 100:
            if (theta != 0.2): 
                snapTheta.append(theta)
                f_d_list.append(f_d)
        theta += (K1 + 2*K2 + 2*K3 + K4)/6
        omega += (L1 + 2*L2 + 2*L3 + L4)/6
        if (theta < -np.pi):                      
            theta += 2*np.pi
        if (theta > np.pi):
            theta -= 2*np.pi
    return snapTheta, f_d_list




def vary_Fd(theta, omega):
    fstep = 0.0005  # We want this to be small
    domain = 0.25   # Set to 0.25 for range 1.35-->1.6
    f_d = 1.35      # Set to 1.35 to start
    size = int(domain/fstep)
    pl.xlabel("Strength of forcing term")
    pl.ylabel("Theta, radians")
    for i in range(size):
        sTheta, fList = chaos_runge_kutta(theta, omega, f_d)
        pl.plot(fList, sTheta, "k.")
        f_d += fstep
    pl.ylim([-np.pi,np.pi])
    pl.show()



# ==================================================================================== #

# Set theta = 0.2 and omega = 0 for initial conditions
vary_Fd(0.2, 0.0) 


