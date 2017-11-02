# -*- coding: utf-8 -*-
"""
Lyapunov Exponent

Created on Fri Oct 27 15:05:11 2017

@author: Maxwell Levin
"""

import numpy as np
import pylab as pl
import random as rd



a = 0       # Start of domain
b = 150     # End of domain
h = 2e-1    # Step size ( keep (b-a)/h within a reasonable size )
g = 9.8     # Acceleration due to gravity
q = 1/2     # Damping constant
m = 1       # Mass of pendulum
L = 9.8     # Length of pendulum
h = 1e-2    # Set our depth/step size
w_d = 2/3   # Set the frequency of our forcing
f_d = 1.2   # Strength of the forcing


tPoints = np.arange(a,b,h)

# forcing function (give it periodic forcing)
def forcing(s):   
    return f_d*np.sin(s)

# Our d[Omega]/dt function
def F2(theta_, w, t):
    return -g*np.sin(theta_)/L - q*w + forcing(w_d*t) 

# Our d[theta]/dt function
def F1(theta_, w, t):
    return w

# Our function that computes the solution to our second order differential
# equation using 4th order Runge-Kutta method
def chaos_runge_kutta(theta, omega=0):
    omegaPoints = []
    thetaPoints = []
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
        theta += (K1 + 2*K2 + 2*K3 + K4)/6
        omega += (L1 + 2*L2 + 2*L3 + L4)/6
        if (theta < -np.pi):                      
            theta += 2*np.pi
        if (theta > np.pi):
            theta -= 2*np.pi
    return thetaPoints


def vary_init_conditions(num_trials):
    log_theta_plot = np.zeros([int((b-a)/h)], float)
    epsi = 1e-6
    for i in range(num_trials):
        theta = rd.uniform(0.2,0.21)
        thetaRand = chaos_runge_kutta(theta, 0)
        thetaEpsi = chaos_runge_kutta(theta + epsi, 0)
        for j in range(len(thetaRand)):
            log_theta_plot[j] += np.log(np.abs(thetaRand[j] - thetaEpsi[j])) + num_trials
            log_theta_plot[j] /= num_trials
    return log_theta_plot


def line_best_fit(arrayX, arrayY):
    # Calculate the mean of x and y
    xbar = 0 
    for x in arrayX:
        xbar += x / len(arrayX)
    ybar = 0
    for y in arrayY:
        ybar += y / len(arrayY)
    # Calculate slope using least squares method
    m = xsum = 0
    for i in range( len(arrayX) ):
        m += (arrayX[i] - xbar)*(arrayY[i] - ybar)
        xsum += (arrayX[i] - xbar)**2
    m /= xsum
    # Plot the line of best fit
    b = ybar - m*xbar
    regression_line = []
    for x in arrayX:
        regression_line.append(m*x+b)
    return regression_line, m
    
    
    


# ======================================================================= #


logPlotA = chaos_runge_kutta(0.2,0)
logPlotB = chaos_runge_kutta(0.2+1e-6,0)

for i in range(len(logPlotA)):
    logPlotA[i] = np.log(np.abs( logPlotA[i] - logPlotB[i] ))

pl.plot(tPoints, logPlotA)
pl.show()

logPlot = vary_init_conditions(100)
pl.plot(tPoints, logPlot, label="Log Plot", color='c')


lineS = 0.05    # Percentage into data to start best fit line
lineF = 0.50    # Percentage into data to finish best fit line
active = int(len(tPoints)*lineS), int(len(tPoints)*lineF)
best_fit, m = line_best_fit(tPoints[active[0]:active[1]], logPlot[active[0]:active[1]])
pl.plot(tPoints[active[0]:active[1]], best_fit, label="Line of best fit", color='r')
pl.legend(loc='lower right')
pl.show()




print("The Lyapunov Exponent is approximately", m, "for this system.\n")