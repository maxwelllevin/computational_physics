# -*- coding: utf-8 -*-
"""
Exercise 8.15 The Double Pendulum

Created on Tue Nov  7 20:53:46 2017

@author: Maxwell
"""


import numpy as np
import pylab as pl


# DONE: Derive an expression for E = T + V in terms of th1, th2, w1, w2
# DONE: Use 4th-Order RK method to solve the equations of motion for L = 40cm,
#       th1=th2=90', w1=w2=0. m1=m2=1kg, graph of E from t=0 to t=100
# TODO: Copy the program and create an animation of the actual pendulums


"""
We seek the solution to the system of first order ODE's:
dth1/dt = w1
dth2/dt = w2
dw1/dt = -w1^2sin(2th1-2th2)+2w2^2sin(th1-th2)+(g/L)[sin(th1-2th2)+3sin(th1)]
                                3-cos(2th1 - 2th2)
dw2/dt =  4w1^2sin(th1-th2)+w2^2sin(2th1-2th2)+2(g/L)[sin(2th1-th2)-sin(th2)]
                                3-cos(2th1 - 2th2)
Which models the motion of a frictionless double pendulum.
Also we seek to solve
E = -mgL(2cos(th1)+costh2) + mL^2(w1^2+0.5w2^2 + w1w2cos(th1-th2))
which should remain constant over time due to nature of our setup. In our
computation we seek to limit dE/dt to below 10^-5 J.
"""


# Our Constants
m = 1
g = 9.8
L = 0.4
h = 1e-3

# Initial Conditions
ith1 = ith2 = np.pi/2
iw1 = iw2 = 0

# Our Domain
tPts = np.arange(0, 100, h)


# Equation of motion -- derivative of theta1
def dth1(th1, th2, w1, w2, time):
    return w1


# Equation of motion -- derivative of theta2
def dth2(th1, th2, w1, w2, time):
    return w2


# Equation of motion -- derivative of omega1
def dw1(th1, th2, w1, w2, time):
    temp = -(w1**2)*np.sin(2*th1-2*th2) + 2*(w1**2)*np.sin(th1-th2)
    temp -= (g/L)*(np.sin(th1-2*th2)+3*np.sin(th1))
    temp /= 3-np.cos(2*th1-2*th2)
    return temp


# Equation of motion -- derivative of omega2
def dw2(th1, th2, w1, w2, time):
    temp = 4*(w1**2)*np.sin(th1-th2) + (w2**2)*np.sin(2*th1-2*th2)
    temp += (2*g/L)*(np.sin(2*th1-th2)-np.sin(th2))
    temp /= 3-np.cos(2*th1-2*th2)
    return temp


# Our function that computes the solution to the second order differential
# equation using 4th order Runge-Kutta method
def doublePendulum(th1, th2, w1, w2):
    # Initial Conditions
    th1Pts = []
    th2Pts = []
    w1Pts = []
    w2Pts = []
    ePts = []

    # RK-step
    for t in tPts:
        th1Pts.append(th1)
        th2Pts.append(th2)
        w1Pts.append(w1)
        w2Pts.append(w2)
        sysEnergy = -m*g*L*(2*np.cos(th1)+np.cos(th2))
        sysEnergy += m*L*L*(w1*w1+0.5*w2*w2+w1*w2*np.cos(th1-th2))
        ePts.append(sysEnergy)

        k1 = h*dth1(th1, th2, w1, w2, t)
        m1 = h*dth2(th1, th2, w1, w2, t)
        n1 = h*dw1(th1, th2, w1, w2, t)
        p1 = h*dw2(th1, th2, w1, w2, t)

        k2 = h*dth1(th1+0.5*k1, th2+0.5*m1, w1+0.5*n1, w2+0.5*p1, t+0.5*h)
        m2 = h*dth2(th1+0.5*k1, th2+0.5*m1, w1+0.5*n1, w2+0.5*p1, t+0.5*h)
        n2 = h*dw1(th1+0.5*k1, th2+0.5*m1, w1+0.5*n1, w2+0.5*p1, t+0.5*h)
        p2 = h*dw2(th1+0.5*k1, th2+0.5*m1, w1+0.5*n1, w2+0.5*p1, t+0.5*h)

        k3 = h*dth1(th1+0.5*k2, th2+0.5*m2, w1+0.5*n2, w2+0.5*p2, t+0.5*h)
        m3 = h*dth2(th1+0.5*k2, th2+0.5*m2, w1+0.5*n2, w2+0.5*p2, t+0.5*h)
        n3 = h*dw1(th1+0.5*k2, th2+0.5*m2, w1+0.5*n2, w2+0.5*p2, t+0.5*h)
        p3 = h*dw2(th1+0.5*k2, th2+0.5*m2, w1+0.5*n2, w2+0.5*p2, t+0.5*h)

        k4 = h*dth1(th1+k3, th2+m3, w1+n3, w2+p3, t+h)
        m4 = h*dth2(th1+k3, th2+m3, w1+n3, w2+p3, t+h)
        n4 = h*dw1(th1+k3, th2+m3, w1+n3, w2+p3, t+h)
        p4 = h*dw2(th1+k3, th2+m3, w1+n3, w2+p3, t+h)

        th1 += (k1 + 2*k2 + 2*k3 + k4)/6
        th2 += (m1 + 2*m2 + 2*m3 + m4)/6
        w1 += (n1 + 2*n2 + 2*n3 + n4)/6
        w2 += (p1 + 2*p2 + 2*p3 + p4)/6

        if (th1 < -np.pi):
            th1 += 2*np.pi
        elif (th1 > np.pi):
            th1 -= 2*np.pi
        if (th2 < -np.pi):
            th2 += 2*np.pi
        elif (th2 > np.pi):
            th2 -= 2*np.pi

    return th1Pts, th2Pts, w1Pts, w2Pts, ePts


# ====================== Main Code ======================== #

th1P, th2P, w1P, w2P, E = doublePendulum(ith1, ith2, iw1, iw2)
#pl.plot(w1P, th1P)
#pl.title("Theta1 vs. Omega1")
#pl.xlabel("Omega1")
#pl.ylabel("Theta1")
#pl.show()


pl.plot(tPts, E)
pl.show()
